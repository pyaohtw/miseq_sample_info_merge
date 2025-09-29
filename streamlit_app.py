# app.py
import streamlit as st
import pandas as pd
from io import BytesIO, StringIO
import datetime
import pytz
import re

# --- CONFIGURATION ---
SHEET_NAME = "Data"
USER_ID_SHEET = "UserID" # === Define UserID sheet name
CORE_COLS = ["Sample_ID", "I7_Index_ID", "index", "I5_Index_ID", "index2", "Amplicon"]
BAD_VALUES = {"#N/A", ""}

# --- PAGE SETUP ---
st.set_page_config(page_title="MiSeq Excel Merger", layout="wide")
st.title("🔗 MiSeq Excel Merger")

# --- DEFAULT PREFIX ---
pst = pytz.timezone("America/Los_Angeles")
today = datetime.datetime.now(pst).date()
default_prefix = f"run{today.year % 100:02d}{today.month:02d}{today.day:02d}mi"
prefix = st.text_input("Filename prefix (defaults to today’s run timestamp)", value=default_prefix)
out_excel = f"{prefix}_sample_info.xlsx"
out_csv = f"{prefix}_miseq.csv"

# --- SESSION STATE ---
state = st.session_state
if "upload_key" not in state:
    state.upload_key = 0
for key in ("expanded_df","raw_df","log_rows","file_combos","cross_dup_combos",
            # === keep per-merge gRNA QC tallies ===
            "grna_qc_logs", "grna_qc_totals",
            # === User ID storage ===
            "user_ids"):
    if key not in state:
        state[key] = None

# --- MERGE & CLEAR BUTTONS ---
col1, col2, _ = st.columns([1,1,4])
with col1:
    merge_clicked = st.button("▶️ Merge")
with col2:
    if st.button("🗑️ Clear uploads"):
        state.upload_key += 1
        for key in ("expanded_df","raw_df","log_rows","file_combos","cross_dup_combos",
                    "grna_qc_logs","grna_qc_totals", "user_ids"):
            state[key] = None

st.markdown("---")

# --- FILE UPLOADER ---
uploader_key = f"uploads_{state.upload_key}"
uploaded_files = st.file_uploader("Upload one or more .xlsx files", type=["xlsx"], accept_multiple_files=True, key=uploader_key)

# === gRNA normalization & QC helpers (unchanged) ===
_GRNA_VALID_RE = re.compile(r"^[ATCG]*$")

def _normalize_grna_val(val):
    """
    Convert a single gRNA value to uppercase string, replace U->T,
    preserve NaN/blank as blank string for QC consistency.
    """
    if pd.isna(val):
        return ""
    s = str(val).strip()
    if s == "":
        return ""
    s = s.upper().replace("U", "T")
    return s

def clean_and_qc_grna(df: pd.DataFrame):
    """
    Returns (df_cleaned, qc_dict)
    - df_cleaned: with gRNA uppercased and U->T applied in-place
    - qc_dict: counts for U->T conversions and remaining invalid gRNA entries
    """
    if "gRNA" not in df.columns:
        # No gRNA column — skip gracefully
        return df, {"u_to_t": 0, "invalid_after": 0, "invalid_examples": []}

    grna_orig = df["gRNA"].copy()

    # Count values containing U/u before normalization (non-empty only)
    contains_u_mask = grna_orig.astype(str).str.contains(r"[Uu]", regex=True, na=False)
    # Normalize
    df["gRNA"] = grna_orig.apply(_normalize_grna_val)

    # Count invalid AFTER normalization (i.e., not strictly A/T/C/G)
    nonempty_mask = df["gRNA"].astype(str).str.len() > 0
    invalid_mask = nonempty_mask & ~df["gRNA"].astype(str).str.match(_GRNA_VALID_RE)
    invalid_after_count = int(invalid_mask.sum())
    invalid_examples = df.loc[invalid_mask, "gRNA"].astype(str).unique().tolist()[:5]

    qc = {
        "u_to_t": int(contains_u_mask.sum()),
        "invalid_after": invalid_after_count,
        "invalid_examples": invalid_examples,
    }
    return df, qc

# --- HELPER: EXPAND BLANK gRNA (unchanged) ---
def expand_gRNA(df):
    non_blank = df[df['gRNA'].notna() & (df['gRNA'].astype(str).str.strip() != '')].copy()
    blank = df[df['gRNA'].isna() | (df['gRNA'].astype(str).str.strip() == '')].copy()
    extras = []
    for _, row in blank.iterrows():
        amp = str(row['Amplicon']).strip()
        matches = non_blank[non_blank['Amplicon'].astype(str).str.strip() == amp]
        for g in matches['gRNA'].dropna().unique():
            new_row = row.copy()
            new_row['gRNA'] = g
            extras.append(new_row)
    df_extra = pd.DataFrame(extras, columns=df.columns) if extras else pd.DataFrame(columns=df.columns)
    return pd.concat([non_blank, df_extra], ignore_index=True), len(blank), len(extras)

# --- PERFORM MERGE & QC ---
if merge_clicked:
    if not uploaded_files:
        st.warning("Please upload at least one file before merging.")
    else:
        raw_list, exp_list, logs = [], [], []
        file_combos = {}
        grna_qc_logs = []
        total_u_to_t = 0
        total_invalid = 0
        total_invalid_examples = set()
        
        user_ids = set() 

        for f in uploaded_files:
            
            df_data = None
            df_user_id = None
            
            # 1. Attempt to read Data sheet (REQUIRED)
            try:
                df_data = pd.read_excel(f, sheet_name=SHEET_NAME, engine='openpyxl')
            except ValueError:
                # Data sheet not found, cannot process this file. Log and skip.
                logs.append({
                    'File': f.name,
                    'Sheet Found': False, 
                    'Input Samples': 0,
                    'Blank gRNA': 0,
                    'gRNA Entries Added': 0,
                    'In-file Dup Combos': 0,
                    'Cross-file Dup': False,
                    'gRNA U→T': 0,
                    'Non-ATCG (post U→T)': 0
                })
                continue
                
            # 2. Attempt to read UserID sheet (OPTIONAL, assuming a header is used)
            try:
                # === UPDATED: Use default header=0 (first row is header)
                df_user_id = pd.read_excel(f, sheet_name=USER_ID_SHEET, engine='openpyxl')
            except ValueError:
                # UserID sheet not found, that's fine, df_user_id remains None
                pass
            
            # --- PROCESS USERID SHEET ---
            if df_user_id is not None:
                
                # Normalize column names to check against (allowing for "UserID" or "User ID")
                cols = [str(c).strip().replace('_', ' ') for c in df_user_id.columns]
                
                # === UPDATED: Look for a column named 'UserID' or 'User ID'
                target_col_index = -1
                if 'UserID'.lower() in [c.lower() for c in cols]:
                    target_col_index = [c.lower() for c in cols].index('userid')
                elif 'User ID'.lower() in [c.lower() for c in cols]:
                    target_col_index = [c.lower() for c in cols].index('user id')
                
                if target_col_index != -1:
                    # Get the actual column name from the original DataFrame
                    target_col_name = df_user_id.columns[target_col_index]
                    
                    if not df_user_id.empty and target_col_name in df_user_id.columns:
                        # Extract data using the column name, ensuring the header is skipped
                        current_ids = df_user_id[target_col_name].astype(str).str.strip().dropna()
                        # Filter out blank IDs and merge
                        user_ids.update(set(id for id in current_ids if id != ""))

            # --- PROCESS DATA SHEET (unchanged) ---
            df = df_data 

            # filter rows by CORE_COLS
            mask = df[CORE_COLS].notna().all(axis=1)
            for col in CORE_COLS:
                mask &= ~df[col].astype(str).isin(BAD_VALUES)
            df_filtered = df[mask].copy()

            # === normalize & QC gRNA before expansion ===
            df_filtered, qc = clean_and_qc_grna(df_filtered)
            total_u_to_t += qc["u_to_t"]
            total_invalid += qc["invalid_after"]
            for ex in qc["invalid_examples"]:
                total_invalid_examples.add(ex)
            grna_qc_logs.append({
                "File": f.name,
                "gRNA U→T": qc["u_to_t"],
                "Non-ATCG (post U→T)": qc["invalid_after"],
                "Examples (up to 5)": ", ".join(qc["invalid_examples"]) if qc["invalid_examples"] else ""
            })

            raw_list.append(df_filtered)

            # track index combos per file
            combos = set(tuple(x) for x in df_filtered[['index','index2']].dropna().apply(tuple, axis=1))
            file_combos[f.name] = combos

            # in-file duplicate combos count
            dup_counts = df_filtered.groupby(['index','index2']).size()
            in_dup_count = sum(1 for c in dup_counts.values if c > 1)

            # expand for Excel
            df_expanded, blank_count, filled_count = expand_gRNA(df_filtered)
            exp_list.append(df_expanded)

            # add log entry
            logs.append({
                'File': f.name,
                'Sheet Found': True, 
                'Input Samples': len(df_filtered),
                'Blank gRNA': blank_count,
                'gRNA Entries Added': filled_count,
                'In-file Dup Combos': in_dup_count,
                'Cross-file Dup': False, 
                'gRNA U→T': qc["u_to_t"],
                'Non-ATCG (post U→T)': qc["invalid_after"]
            })

        # ... (rest of log and state storage logic remains the same) ...
        # cross-file duplicate detection
        combo_files_map = {}
        for fname, combos in file_combos.items():
            for combo in combos:
                combo_files_map.setdefault(combo, []).append(fname)
        cross_dup = {c: fs for c,fs in combo_files_map.items() if len(fs)>1}

        # mark cross-file dup
        for entry in logs:
            entry['Cross-file Dup'] = any(entry['File'] in fs for fs in cross_dup.values())

        # store data and logs
        state.raw_df = pd.concat(raw_list, ignore_index=True) if raw_list else pd.DataFrame()
        state.expanded_df = pd.concat(exp_list, ignore_index=True) if exp_list else pd.DataFrame()

        log_df = pd.DataFrame(logs)
        # sum numeric cols safely
        totals = {}
        for col in log_df.columns:
            if col in ('File', 'Sheet Found', 'Cross-file Dup'):
                totals[col] = '' if col != 'File' else 'Total'
            else:
                try:
                    totals[col] = int(pd.to_numeric(log_df[col], errors='coerce').fillna(0).sum())
                except Exception:
                    totals[col] = ''
        state.log_rows = pd.concat([log_df, pd.DataFrame([totals])], ignore_index=True)

        state.file_combos = file_combos
        state.cross_dup_combos = cross_dup

        # === save gRNA QC summary for notifications
        state.grna_qc_logs = pd.DataFrame(grna_qc_logs) if grna_qc_logs else None
        state.grna_qc_totals = {
            "u_to_t": total_u_to_t,
            "invalid_after": total_invalid,
            "invalid_examples": sorted(list(total_invalid_examples))[:5]
        }
        
        # === save merged User IDs
        state.user_ids = sorted(list(user_ids))

# --- DISPLAY RESULTS & DOWNLOADS ---
if state.log_rows is not None:
    # === show gRNA QC notifications ===
    if state.grna_qc_totals:
        if state.grna_qc_totals["u_to_t"] > 0:
            st.info(f"gRNA cleanup: Converted **{state.grna_qc_totals['u_to_t']}** entries from U→T (case-insensitive).")
        if state.grna_qc_totals["invalid_after"] > 0:
            examples = ", ".join(state.grna_qc_totals["invalid_examples"])
            msg = f"Found **{state.grna_qc_totals['invalid_after']}** gRNA entr{'y' if state.grna_qc_totals['invalid_after']==1 else 'ies'} containing non-A/T/C/G characters **after** U→T conversion."
            if examples:
                msg += f" Examples: {examples}"
            st.warning(msg)

    st.subheader('Merge Log')
    st.table(state.log_rows)

    # === optional per-file gRNA QC table (collapsed) ===
    if state.grna_qc_logs is not None and not state.grna_qc_logs.empty:
        with st.expander("Show per-file gRNA QC details"):
            st.table(state.grna_qc_logs)
            
    # === show merged User IDs
    if state.user_ids:
        user_id_str = ", ".join(state.user_ids)
        st.info(f"Merged Unique User IDs: {user_id_str}")

    if state.cross_dup_combos:
        combos_str = ', '.join(f"{i}/{j}" for i,j in state.cross_dup_combos)
        st.warning(f"Cross-file duplicate index combos across files: {combos_str}")

    # Excel download
    excel_buf = BytesIO()
    with pd.ExcelWriter(excel_buf, engine='openpyxl') as w:
        state.expanded_df.to_excel(w, index=False, sheet_name=SHEET_NAME)
        
        # === Write UserID sheet
        if state.user_ids:
            # Create a DataFrame with a single column 'UserID' for consistent output
            user_df = pd.DataFrame(
                {
                    'UserID': [", ".join(state.user_ids)]
                }
            )
            # Write to Excel with the 'UserID' header
            user_df.to_excel(w, index=False, sheet_name=USER_ID_SHEET)
            
    excel_buf.seek(0)
    st.download_button('📥 Download merged Excel', data=excel_buf, file_name=out_excel,
                       mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')

    # CSV download
    csv_buf = StringIO()
    csv_buf.write('[Header]\n')
    csv_buf.write(f'Experiment Name,{prefix}\n')
    csv_buf.write(f'Date,{today.month}/{today.day}/{today.year}\n')
    csv_buf.write('Workflow,GenerateFASTQ\n')
    csv_buf.write('[Reads]\n300\n')
    csv_buf.write('[Settings]\n')
    csv_buf.write('[Data]\n')
    csv_buf.write('Sample_ID,Sample_Name,I7_Index_ID,index,I5_INDEX_ID,index2,Sample_Project,Description\n')
    for _,row in state.raw_df.iterrows():
        csv_buf.write(
            f"{row['Sample_ID']},,{row['I7_Index_ID']},{row['index']}"
            f",{row['I5_Index_ID']},{row['index2']},,\n"
        )
    data_c = csv_buf.getvalue().encode('utf-8')
    st.download_button('📥 Download Miseq CSV', data=data_c, file_name=out_csv, mime='text/csv')

    st.subheader('Merged Data Preview')
    st.dataframe(state.expanded_df, use_container_width=True)

else:
    if merge_clicked:
        st.error('No valid rows found.')