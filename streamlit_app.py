# app.py
import streamlit as st
import pandas as pd
from io import BytesIO, StringIO
import datetime

# --- CONFIGURATION ---
SHEET_NAME = "Data"
CORE_COLS = ["Sample_ID", "I7_Index_ID", "index", "I5_Index_ID", "index2", "Amplicon"]
BAD_VALUES = {"#N/A", ""}

# --- PAGE SETUP ---
st.set_page_config(page_title="MiSeq Excel Merger", layout="wide")
st.title("🔗 MiSeq Excel Merger")

# --- DEFAULT PREFIX ---
today = datetime.date.today()
default_prefix = f"run{today.year % 100:02d}{today.month:02d}{today.day:02d}mi"
prefix = st.text_input("Filename prefix (defaults to today’s run timestamp)", value=default_prefix)
out_excel = f"{prefix}_sample_info.xlsx"
out_csv = f"{prefix}_miseq.csv"

# --- SESSION STATE ---
state = st.session_state
if "upload_key" not in state:
    state.upload_key = 0
for key in ("expanded_df","raw_df","log_rows","file_combos","cross_dup_combos"): 
    if key not in state:
        state[key] = None

# --- MERGE & CLEAR BUTTONS ---
col1, col2, _ = st.columns([1,1,4])
with col1:
    merge_clicked = st.button("▶️ Merge")
with col2:
    if st.button("🗑️ Clear uploads"):
        state.upload_key += 1
        for key in ("expanded_df","raw_df","log_rows","file_combos","cross_dup_combos"): 
            state[key] = None

st.markdown("---")

# --- FILE UPLOADER ---
uploader_key = f"uploads_{state.upload_key}"
uploaded_files = st.file_uploader("Upload one or more .xlsx files", type=["xlsx"], accept_multiple_files=True, key=uploader_key)

# --- HELPER: EXPAND BLANK gRNA ---
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
        for f in uploaded_files:
            # attempt to read Data sheet
            try:
                df = pd.read_excel(f, sheet_name=SHEET_NAME, engine='openpyxl')
                sheet_found = True
            except ValueError:
                # log missing sheet and skip
                logs.append({
                    'File': f.name,
                    'Sheet Found': False,
                    'Input Samples': 0,
                    'Blank gRNA': 0,
                    'gRNA Entries Added': 0,
                    'In-file Dup Combos': 0,
                    'Cross-file Dup': False
                })
                continue
            # proceed with QC and merge
            # filter rows by CORE_COLS
            mask = df[CORE_COLS].notna().all(axis=1)
            for col in CORE_COLS:
                mask &= ~df[col].astype(str).isin(BAD_VALUES)
            df_filtered = df[mask].copy()
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
                'Cross-file Dup': False  # annotate later
            })
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
        totals = {col: ('' if col=='File' or col=='Sheet Found' else int(log_df[col].sum())) for col in log_df.columns}
        totals['File'] = 'Total'
        totals['Sheet Found'] = ''
        state.log_rows = pd.concat([log_df, pd.DataFrame([totals])], ignore_index=True)
        state.file_combos = file_combos
        state.cross_dup_combos = cross_dup

# --- DISPLAY RESULTS & DOWNLOADS ---
if state.log_rows is not None:
    st.subheader('Merge Log')
    st.table(state.log_rows)
    if state.cross_dup_combos:
        combos_str = ', '.join(f"{i}/{j}" for i,j in state.cross_dup_combos)
        st.warning(f"Cross-file duplicate index combos across files: {combos_str}")
    # Excel download
    excel_buf = BytesIO()
    with pd.ExcelWriter(excel_buf, engine='openpyxl') as w:
        state.expanded_df.to_excel(w, index=False, sheet_name=SHEET_NAME)
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
