# app.py
import streamlit as st
import pandas as pd
from io import BytesIO, StringIO
import datetime
import pytz
import re
from openpyxl import load_workbook
from openpyxl.styles import numbers

# --- CONFIGURATION ---
SHEET_NAME = "Data"
USER_ID_SHEET = "UserID"
CORE_COLS = ["Sample_ID", "I7_Index_ID", "index", "I5_Index_ID", "index2", "Amplicon"]
BAD_VALUES = {"#N/A", ""}

# Columns to treat as text (never numeric)
TEXT_BLANK_COLS = ["ELN_ID", "Isoform_Sample_ID", "PAM"]

# Final output columns for CRISPResso
FINAL_COLS = [
    "Sample_ID", "Sample_Name", "I7_Index_ID", "index", "I5_Index_ID", "index2",
    "Sample_Project", "Description", "ELN_ID", "Isoform_Sample_ID", "PAM",
    "gRNA", "Amplicon", "Exon",
    "Expected_HDR_Amplicon", "Quantification_Window_Coordinates",
    "Quantification_Window_Center", "Plot_Window_Size", "ngRNA"
]

NGRNA_RAW_HEADERS = ["ngRNA\n(nicking RNA)", "ngRNA (nicking RNA)"]
NGRNA_FINAL_HEADER = "ngRNA"

# --- PAGE SETUP ---
st.set_page_config(page_title="MiSeq Excel Merger", layout="wide")
st.title("🔗 MiSeq Excel Merger")

# --- DEFAULT PREFIX ---
pst = pytz.timezone("America/Los_Angeles")
today = datetime.datetime.now(pst).date()
default_prefix = f"run{today.year % 100:02d}{today.month:02d}{today.day:02d}mi"
prefix = st.text_input("Filename prefix (defaults to today's run timestamp)", value=default_prefix)
out_excel = f"{prefix}_sample_info.xlsx"
out_csv = f"{prefix}_miseq.csv"

# --- SESSION STATE ---
state = st.session_state
if "upload_key" not in state:
    state.upload_key = 0
for key in ("expanded_df", "raw_df", "log_rows", "file_combos", "cross_dup_combos",
            "grna_qc_logs", "grna_qc_totals", "user_ids"):
    if key not in state:
        state[key] = None

# --- MERGE & CLEAR BUTTONS ---
col1, col2, _ = st.columns([1, 1, 4])
with col1:
    merge_clicked = st.button("▶️ Merge")
with col2:
    if st.button("🗑️ Clear uploads"):
        state.upload_key += 1
        for key in ("expanded_df", "raw_df", "log_rows", "file_combos", "cross_dup_combos",
                    "grna_qc_logs", "grna_qc_totals", "user_ids"):
            state[key] = None

st.markdown("---")

# --- FILE UPLOADER ---
uploaded_files = st.file_uploader(
    "Upload one or more .xlsx files",
    type=["xlsx"],
    accept_multiple_files=True,
    key=f"uploads_{state.upload_key}"
)

# --- HELPER FUNCTIONS ---
_GRNA_VALID_RE = re.compile(r"^[ATCG]*$")

def normalize_grna_val(val):
    """Convert gRNA to uppercase, replace U->T, preserve blank as empty string."""
    if pd.isna(val):
        return ""
    s = str(val).strip()
    if s == "":
        return ""
    return s.upper().replace("U", "T")

def clean_cell(x):
    """Convert NaN/None/null strings to empty string."""
    if x is None or pd.isna(x):
        return ""
    s = str(x).strip()
    if s.lower() in ("nan", "none", "null", "0"):
        return ""
    return s

def blankify(df):
    """Enforce true blanks everywhere - convert NaN and 'nan' strings to empty strings."""
    if df is None or df.empty:
        return df
    df = df.where(pd.notna(df), "")
    return df.applymap(clean_cell)

def clean_and_qc_grna(df):
    """Returns (df_cleaned, qc_dict) with gRNA normalized and QC stats."""
    if "gRNA" not in df.columns:
        return df, {"u_to_t": 0, "invalid_after": 0, "invalid_examples": []}

    grna_orig = df["gRNA"].copy()
    contains_u_mask = grna_orig.astype(str).str.contains(r"[Uu]", regex=True, na=False)
    
    df["gRNA"] = grna_orig.apply(normalize_grna_val)
    
    nonempty_mask = df["gRNA"].astype(str).str.len() > 0
    invalid_mask = nonempty_mask & ~df["gRNA"].astype(str).str.match(_GRNA_VALID_RE)
    invalid_after_count = int(invalid_mask.sum())
    invalid_examples = df.loc[invalid_mask, "gRNA"].astype(str).unique().tolist()[:5]

    return df, {
        "u_to_t": int(contains_u_mask.sum()),
        "invalid_after": invalid_after_count,
        "invalid_examples": invalid_examples,
    }

def expand_gRNA(df):
    """Expand blank gRNA rows by copying non-blank gRNA values with same Amplicon."""
    non_blank = df[df["gRNA"].notna() & (df["gRNA"].astype(str).str.strip() != "")].copy()
    blank = df[df["gRNA"].isna() | (df["gRNA"].astype(str).str.strip() == "")].copy()
    
    extras = []
    for _, row in blank.iterrows():
        amp = str(row["Amplicon"]).strip()
        matches = non_blank[non_blank["Amplicon"].astype(str).str.strip() == amp]
        for g in matches["gRNA"].dropna().unique():
            new_row = row.copy()
            new_row["gRNA"] = g
            extras.append(new_row)
    
    df_extra = pd.DataFrame(extras, columns=df.columns) if extras else pd.DataFrame(columns=df.columns)
    return pd.concat([non_blank, df_extra], ignore_index=True), len(blank), len(extras)

def ensure_columns_and_order(df, columns):
    """Ensure all required columns exist and are in correct order."""
    df_out = df.copy()
    
    for col in columns:
        if col not in df_out.columns:
            df_out[col] = ""
    
    df_out = df_out[columns]
    df_out = blankify(df_out)
    
    # Extra guard for text columns - ensure they stay as strings
    for c in TEXT_BLANK_COLS:
        if c in df_out.columns:
            df_out[c] = df_out[c].apply(clean_cell)
    
    return df_out

def write_excel_with_blanks(df, excel_buf, sheet_name, user_ids=None):
    """Write Excel ensuring truly blank cells (not 0) for empty values."""
    # Replace empty strings with None for pandas
    df_to_write = df.replace("", None)
    
    # Write initial Excel
    with pd.ExcelWriter(excel_buf, engine="openpyxl") as writer:
        df_to_write.to_excel(writer, index=False, sheet_name=sheet_name)
        if user_ids:
            user_df = pd.DataFrame({"UserID": [", ".join(user_ids)]})
            user_df.to_excel(writer, index=False, sheet_name=USER_ID_SHEET)
    
    # Reopen to format cells as text (prevents 0 display)
    excel_buf.seek(0)
    wb = load_workbook(excel_buf)
    ws = wb[sheet_name]
    
    # Set all data cells to text format
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row, min_col=1, max_col=ws.max_column):
        for cell in row:
            cell.number_format = numbers.FORMAT_TEXT
            # If cell has value 0 or None, clear it
            if cell.value == 0 or cell.value is None:
                cell.value = None
    
    # Save back to buffer
    excel_buf.seek(0)
    excel_buf.truncate()
    wb.save(excel_buf)
    excel_buf.seek(0)

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
            # Read Data sheet (REQUIRED)
            try:
                df_data = pd.read_excel(
                    f, sheet_name=SHEET_NAME, engine="openpyxl",
                    dtype=str, keep_default_na=False, na_filter=False
                )
            except ValueError:
                logs.append({
                    "File": f.name, "Sheet Found": False, "Input Samples": 0,
                    "Blank gRNA": 0, "gRNA Entries Added": 0, "In-file Dup Combos": 0,
                    "Cross-file Dup": False, "gRNA U→T": 0, "Non-ATCG (post U→T)": 0
                })
                continue

            # Read UserID sheet (OPTIONAL)
            try:
                df_user_id = pd.read_excel(
                    f, sheet_name=USER_ID_SHEET, engine="openpyxl",
                    dtype=str, keep_default_na=False, na_filter=False
                )
                cols_lower = [str(c).strip().replace("_", " ").lower() for c in df_user_id.columns]
                target_idx = next((i for i, c in enumerate(cols_lower) if c in ["userid", "user id"]), -1)
                
                if target_idx != -1:
                    target_col = df_user_id.columns[target_idx]
                    current_ids = df_user_id[target_col].astype(str).str.strip()
                    user_ids.update(i for i in current_ids if i and i != "")
            except ValueError:
                pass

            # Process Data sheet
            df = blankify(df_data)
            
            # Normalize ngRNA column header
            for h in NGRNA_RAW_HEADERS:
                if h in df.columns and NGRNA_FINAL_HEADER not in df.columns:
                    df = df.rename(columns={h: NGRNA_FINAL_HEADER})
                    break

            # Filter rows by CORE_COLS
            mask = df[CORE_COLS].notna().all(axis=1)
            for col in CORE_COLS:
                mask &= ~df[col].astype(str).isin(BAD_VALUES)
            df_filtered = blankify(df[mask].copy())

            # Normalize & QC gRNA
            df_filtered, qc = clean_and_qc_grna(df_filtered)
            df_filtered = blankify(df_filtered)

            total_u_to_t += qc["u_to_t"]
            total_invalid += qc["invalid_after"]
            total_invalid_examples.update(qc["invalid_examples"])
            grna_qc_logs.append({
                "File": f.name, "gRNA U→T": qc["u_to_t"],
                "Non-ATCG (post U→T)": qc["invalid_after"],
                "Examples (up to 5)": ", ".join(qc["invalid_examples"]) if qc["invalid_examples"] else ""
            })

            raw_list.append(df_filtered)

            # Track index combos per file
            combos = set(tuple(x) for x in df_filtered[["index", "index2"]].dropna().apply(tuple, axis=1))
            file_combos[f.name] = combos

            # In-file duplicate combos
            dup_counts = df_filtered.groupby(["index", "index2"]).size()
            in_dup_count = sum(1 for c in dup_counts.values if c > 1)

            # Expand for Excel
            df_expanded, blank_count, filled_count = expand_gRNA(df_filtered)
            df_expanded = blankify(df_expanded)
            exp_list.append(df_expanded)

            logs.append({
                "File": f.name, "Sheet Found": True, "Input Samples": len(df_filtered),
                "Blank gRNA": blank_count, "gRNA Entries Added": filled_count,
                "In-file Dup Combos": in_dup_count, "Cross-file Dup": False,
                "gRNA U→T": qc["u_to_t"], "Non-ATCG (post U→T)": qc["invalid_after"]
            })

        # Cross-file duplicate detection
        combo_files_map = {}
        for fname, combos in file_combos.items():
            for combo in combos:
                combo_files_map.setdefault(combo, []).append(fname)
        cross_dup = {c: fs for c, fs in combo_files_map.items() if len(fs) > 1}

        for entry in logs:
            entry["Cross-file Dup"] = any(entry["File"] in fs for fs in cross_dup.values())

        # Store data
        state.raw_df = blankify(pd.concat(raw_list, ignore_index=True)) if raw_list else pd.DataFrame()
        state.expanded_df = blankify(pd.concat(exp_list, ignore_index=True)) if exp_list else pd.DataFrame()

        # Create log with totals
        log_df = pd.DataFrame(logs)
        totals = {"File": "Total"}
        for col in log_df.columns:
            if col not in ("File", "Sheet Found", "Cross-file Dup"):
                try:
                    totals[col] = int(pd.to_numeric(log_df[col], errors="coerce").fillna(0).sum())
                except Exception:
                    totals[col] = ""
        state.log_rows = pd.concat([log_df, pd.DataFrame([totals])], ignore_index=True)

        state.file_combos = file_combos
        state.cross_dup_combos = cross_dup
        state.grna_qc_logs = pd.DataFrame(grna_qc_logs) if grna_qc_logs else None
        state.grna_qc_totals = {
            "u_to_t": total_u_to_t,
            "invalid_after": total_invalid,
            "invalid_examples": sorted(list(total_invalid_examples))[:5]
        }
        state.user_ids = sorted(list(user_ids))

# --- DISPLAY RESULTS & DOWNLOADS ---
if state.log_rows is not None:
    # Show gRNA QC notifications
    if state.grna_qc_totals:
        if state.grna_qc_totals["u_to_t"] > 0:
            st.info(f"gRNA cleanup: Converted **{state.grna_qc_totals['u_to_t']}** entries from U→T (case-insensitive).")
        if state.grna_qc_totals["invalid_after"] > 0:
            examples = ", ".join(state.grna_qc_totals["invalid_examples"])
            msg = f"Found **{state.grna_qc_totals['invalid_after']}** gRNA entry(ies) with non-A/T/C/G characters after U→T conversion."
            if examples:
                msg += f" Examples: {examples}"
            st.warning(msg)

    st.subheader("Merge Log")
    st.table(state.log_rows)

    if state.grna_qc_logs is not None and not state.grna_qc_logs.empty:
        with st.expander("Show per-file gRNA QC details"):
            st.table(state.grna_qc_logs)

    if state.user_ids:
        st.info(f"Merged Unique User IDs: {', '.join(state.user_ids)}")

    if state.cross_dup_combos:
        combos_str = ", ".join(f"{i}/{j}" for i, j in state.cross_dup_combos)
        st.warning(f"Cross-file duplicate index combos: {combos_str}")

    # Prepare final dataframe with proper column order
    final_df = ensure_columns_and_order(state.expanded_df, FINAL_COLS)

    # Warning for originally blank gRNA rows
    if "gRNA" in state.raw_df.columns:
        orig_blank_mask = state.raw_df["gRNA"].astype(str).str.strip().eq("")
        if orig_blank_mask.any():
            n_blank = int(orig_blank_mask.sum())
            id_col = "Sample_ID" if "Sample_ID" in state.raw_df.columns else None
            examples = (
                state.raw_df.loc[orig_blank_mask, id_col].astype(str).head(5).tolist()
                if id_col else []
            )
            if examples:
                st.warning(f"⚠️ {n_blank} row(s) had empty gRNA and were auto-filled. Examples: {', '.join(examples)}")

    # Warning for Expected_HDR_Amplicon without Quantification_Window_Coordinates
    hdr_mask = final_df["Expected_HDR_Amplicon"].astype(str).str.strip().ne("")
    qwc_blank = final_df["Quantification_Window_Coordinates"].astype(str).str.strip().eq("")
    prob_mask = hdr_mask & qwc_blank
    if prob_mask.any():
        n_prob = int(prob_mask.sum())
        examples = final_df.loc[prob_mask, "Sample_ID"].astype(str).head(5).tolist()
        st.warning(f"⚠️ {n_prob} row(s) have Expected_HDR_Amplicon but missing Quantification_Window_Coordinates. Examples: {', '.join(examples)}")

    # Excel download with proper blank handling
    excel_buf = BytesIO()
    write_excel_with_blanks(final_df, excel_buf, SHEET_NAME, state.user_ids)
    
    st.download_button(
        "📥 Download merged Excel (CRISPResso)",
        data=excel_buf,
        file_name=out_excel,
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

    # CSV download
    csv_buf = StringIO()
    csv_buf.write("[Header]\n")
    csv_buf.write(f"Experiment Name,{prefix}\n")
    csv_buf.write(f"Date,{today.month}/{today.day}/{today.year}\n")
    csv_buf.write("Workflow,GenerateFASTQ\n")
    csv_buf.write("[Reads]\n300\n")
    csv_buf.write("[Settings]\n")
    csv_buf.write("[Data]\n")
    csv_buf.write("Sample_ID,Sample_Name,I7_Index_ID,index,I5_INDEX_ID,index2,Sample_Project,Description\n")
    for _, row in state.raw_df.iterrows():
        csv_buf.write(
            f"{row['Sample_ID']},,{row['I7_Index_ID']},{row['index']}"
            f",{row['I5_Index_ID']},{row['index2']},,\n"
        )
    st.download_button(
        "📥 Download Miseq CSV",
        data=csv_buf.getvalue().encode("utf-8"),
        file_name=out_csv,
        mime="text/csv"
    )

    st.subheader("Merged Data Preview")
    st.dataframe(state.expanded_df, use_container_width=True)

else:
    if merge_clicked:
        st.error("No valid rows found.")