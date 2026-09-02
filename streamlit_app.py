import datetime
import re
from io import BytesIO, StringIO

import pandas as pd
import pytz
import streamlit as st
from openpyxl import load_workbook
from openpyxl.styles import numbers


SHEET_NAME = "Data"
USER_ID_SHEET = "UserID"
CORE_COLS = ["Sample_ID", "I7_Index_ID", "index", "I5_Index_ID", "index2", "Amplicon", "gRNA"]
MISEQ_REQUIRED_COLS = ["Sample_ID", "I7_Index_ID", "index", "I5_Index_ID", "index2"]
TEXT_BLANK_COLS = ["ELN_ID", "Isoform_Sample_ID", "PAM", "Base_Editing_Type", "BE_Q30_cutoff"]
FINAL_COLS = [
    "Sample_ID", "Sample_Name", "I7_Index_ID", "index", "I5_Index_ID", "index2",
    "Sample_Project", "Description", "ELN_ID", "Isoform_Sample_ID", "PAM",
    "gRNA", "Amplicon", "Exon",
    "Expected_HDR_Amplicon", "Quantification_Window_Coordinates",
    "Quantification_Window_Center", "Plot_Window_Size", "ngRNA", "Base_Editing_Type", "BE_Q30_cutoff"
]
NGRNA_RAW_HEADERS = ["ngRNA\n(nicking RNA)", "ngRNA (nicking RNA)"]
NGRNA_FINAL_HEADER = "ngRNA"
BASE_EDITING_COL = "Base_Editing_Type"
BE_Q30_CUTOFF_COL = "BE_Q30_cutoff"
VALID_BASE_EDITING_TYPES = {"ABE", "CBE", "BOTH"}
VALID_BE_Q30_CUTOFF_VALUES = {"ON", "OFF"}
BASE_EDITING_CONFLICT_COLS = [
    "Expected_HDR_Amplicon",
    "Quantification_Window_Coordinates",
    "Quantification_Window_Center",
    "Plot_Window_Size",
    NGRNA_FINAL_HEADER,
]
# CRISPResso excludes the terminal 15 bp from the right side of the
# quantification window.  The merge output therefore caps only the coordinate
# end at the last eligible reference position; the coordinate start remains
# allowed anywhere within the retained amplicon.
CRISPRESSO_EXCLUDE_BP_FROM_RIGHT = 15
# Keep a 2-bp safety margin when a symmetric plot window reaches either
# boundary of the trimmed reference.  CRISPResso can reject an exact-boundary
# case such as cut point 160 + plot window 55 == reference length 215.
PLOT_WINDOW_SAFETY_MARGIN = 2
# Defaults used only to assess rows for base-editing/cutting samples whose
# quantification-window fields are intentionally blank.
BLANK_WINDOW_QUANTIFICATION_CENTER = -10
BLANK_WINDOW_PLOT_SIZE = 22
# User-requested CRISPResso right-side exclusion used by this preflight check.
BLANK_WINDOW_CRISPRESSO_EXCLUSION = 5
DNA_ONLY_RE = re.compile(r"^[ATCG]*$")
DNA_OR_U_RE = re.compile(r"^[ATCGU]*$")
WELL_SUFFIX_RE = re.compile(r"([A-H](?:[1-9]|1[0-2]))$", re.IGNORECASE)
QC_OPTIONAL_COLS = [
    "Exon", "Expected_HDR_Amplicon", "Quantification_Window_Coordinates",
    "Quantification_Window_Center", "Plot_Window_Size", BASE_EDITING_COL, NGRNA_FINAL_HEADER,
]
QC_ANCHOR_COLS = ["Sample_ID", "gRNA", "Amplicon", "Exon", "Expected_HDR_Amplicon", BASE_EDITING_COL]
SEVERITY_LABELS = {
    "Error": "Error, must fix❗",
    "Warning": "Warning, review⚠️",
    "Info": "Info, awarenessℹ️",
}
SEVERITY_ORDER = {"Error": 0, "Warning": 1, "Info": 2}


st.set_page_config(page_title="CRISPResso QC / MiSeq Merge Tool", layout="wide")
st.title("🧪 CRISPResso QC / MiSeq Merge Tool")

pst = pytz.timezone("America/Los_Angeles")
today = datetime.datetime.now(pst).date()
default_prefix = f"run{today.year % 100:02d}{today.month:02d}{today.day:02d}mi"

mode = st.radio(
    "Mode",
    options=["QC Screening", "File Merge"],
    index=0,
    horizontal=True,
    help="QC Screening is read-only and reports issues. File Merge generates CRISPResso and MiSeq outputs.",
)

prefix = st.text_input("Filename prefix (used for File Merge outputs)", value=default_prefix)
out_excel = f"{prefix}_sample_info.xlsx"
out_csv = f"{prefix}_miseq.csv"

show_results_as_table = st.checkbox(
    "Display results as interactive table (unchecked = code/text view)",
    value=False,
    key="display_as_table",
    help=(
        "Master switch for all result previews. When on, every table below renders as an "
        "interactive dataframe; when off (default), they render as plain code/text blocks."
    ),
)

state = st.session_state
if "upload_key" not in state:
    state.upload_key = 0
for key in ("qc_results", "merge_results"):
    if key not in state:
        state[key] = None
if "miseq_reads" not in state:
    state.miseq_reads = 284
if "miseq_reads_slider" not in state:
    state.miseq_reads_slider = 284
if "miseq_reads_input" not in state:
    state.miseq_reads_input = 284


def sync_miseq_reads_from_slider():
    value = int(st.session_state.miseq_reads_slider_widget)
    st.session_state.miseq_reads = value
    st.session_state.miseq_reads_slider = value
    st.session_state.miseq_reads_input = value
    st.session_state.miseq_reads_input_widget = value



def sync_miseq_reads_from_input():
    value = int(st.session_state.miseq_reads_input_widget)
    st.session_state.miseq_reads = value
    st.session_state.miseq_reads_input = value
    st.session_state.miseq_reads_slider = value
    st.session_state.miseq_reads_slider_widget = value



def reset_miseq_reads():
    st.session_state.miseq_reads = 284
    st.session_state.miseq_reads_slider = 284
    st.session_state.miseq_reads_input = 284
    st.session_state.miseq_reads_slider_widget = 284
    st.session_state.miseq_reads_input_widget = 284


auto_fix_sample_id_dups = False
if mode == "File Merge":
    auto_fix_sample_id_dups = st.checkbox(
        "Auto-correct duplicate cleaned Sample_IDs",
        value=True,
        help=(
            "When enabled, duplicate Sample_ID values (after removing non-alphanumeric characters) are "
            "detected across ALL uploaded files pooled together, then auto-corrected. The first occurrence "
            "keeps its ID; each later duplicate gets an underscore and occurrence number appended at the "
            "end (e.g. AD1G12, AD1G12_2, AD1G12_3)."
        ),
    )
    slider_col, input_col, reset_col = st.columns([3, 1.2, 1])
    with slider_col:
        if "miseq_reads_slider_widget" not in st.session_state:
            st.session_state.miseq_reads_slider_widget = int(st.session_state.get("miseq_reads", 284))
        st.slider(
            "MiSeq read length (single-end)",
            min_value=25,
            max_value=300,
            key="miseq_reads_slider_widget",
            on_change=sync_miseq_reads_from_slider,
            help="Value written under [Reads] in the MiSeq CSV.",
        )
    with input_col:
        if "miseq_reads_input_widget" not in st.session_state:
            st.session_state.miseq_reads_input_widget = int(st.session_state.get("miseq_reads", 284))
        st.number_input(
            "Manual input",
            min_value=25,
            max_value=300,
            step=1,
            key="miseq_reads_input_widget",
            on_change=sync_miseq_reads_from_input,
        )
    with reset_col:
        st.write("")
        st.write("")
        st.button("Set to 284", on_click=reset_miseq_reads)

col1, col2, _ = st.columns([1, 1, 4])
with col1:
    run_clicked = st.button("▶️ Run QC" if mode == "QC Screening" else "▶️ Merge")
with col2:
    if st.button("🗑️ Clear uploads"):
        state.upload_key += 1
        state.qc_results = None
        state.merge_results = None

st.markdown("---")

uploaded_files = st.file_uploader(
    "Upload one or more .xlsx files",
    type=["xlsx"],
    accept_multiple_files=True,
    key=f"uploads_{state.upload_key}",
)


def normalize_header_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for raw_name in NGRNA_RAW_HEADERS:
        if raw_name in df.columns and NGRNA_FINAL_HEADER not in df.columns:
            df = df.rename(columns={raw_name: NGRNA_FINAL_HEADER})
    return df


def clean_cell_minimal(value) -> str:
    if value is None or pd.isna(value):
        return ""
    s = str(value).strip()
    if s.lower() in {"nan", "none", "null", "0", "0.0"}:
        return ""
    return s


def blankify_minimal(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return df
    df = df.where(pd.notna(df), "")
    return df.map(clean_cell_minimal)


def clean_cell_merge(value) -> str:
    if value is None or pd.isna(value):
        return ""
    s = str(value).strip()
    if s.lower() in {"nan", "none", "null", "0", "0.0"}:
        return ""
    return s


def blankify_merge(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return df
    df = df.where(pd.notna(df), "")
    return df.map(clean_cell_merge)


def normalize_dna(value) -> str:
    return clean_cell_minimal(value).upper()


def normalize_grna_for_merge(value) -> str:
    s = clean_cell_merge(value)
    if not s:
        return ""
    return s.upper().replace("U", "T")


def normalize_base_editing_type(value) -> str:
    s = clean_cell_minimal(value)
    return s.upper() if s else ""


def normalize_be_q30_cutoff(value) -> str:
    s = clean_cell_minimal(value)
    return s.upper() if s else ""


def reverse_complement(seq: str) -> str:
    s = clean_cell_minimal(seq).upper().replace("U", "T")
    trans = str.maketrans("ATCG", "TAGC")
    return s.translate(trans)[::-1]


def ensure_columns_and_order(df: pd.DataFrame, columns) -> pd.DataFrame:
    df = df.copy()
    for col in columns:
        if col not in df.columns:
            df[col] = ""
    df = df[columns]
    df = blankify_merge(df)
    for col in TEXT_BLANK_COLS:
        if col in df.columns:
            df[col] = df[col].apply(clean_cell_merge)
    return df


def write_excel_with_blanks(df: pd.DataFrame, excel_buf: BytesIO, sheet_name: str, user_ids=None):
    df_to_write = df.replace("", None)
    with pd.ExcelWriter(excel_buf, engine="openpyxl") as writer:
        df_to_write.to_excel(writer, index=False, sheet_name=sheet_name)
        if user_ids:
            pd.DataFrame({"UserID": [", ".join(user_ids)]}).to_excel(writer, index=False, sheet_name=USER_ID_SHEET)

    excel_buf.seek(0)
    wb = load_workbook(excel_buf)
    ws = wb[sheet_name]
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row, min_col=1, max_col=ws.max_column):
        for cell in row:
            cell.number_format = numbers.FORMAT_TEXT
            if cell.value == 0 or cell.value is None:
                cell.value = None

    excel_buf.seek(0)
    excel_buf.truncate()
    wb.save(excel_buf)
    excel_buf.seek(0)


def read_user_ids(excel_file) -> list:
    user_ids = set()
    try:
        df_user_id = pd.read_excel(
            excel_file,
            sheet_name=USER_ID_SHEET,
            engine="openpyxl",
            dtype=str,
            keep_default_na=False,
            na_filter=False,
        )
    except ValueError:
        return []

    cols_lower = [str(c).strip().replace("_", " ").lower() for c in df_user_id.columns]
    target_idx = next((i for i, c in enumerate(cols_lower) if c in {"userid", "user id"}), -1)
    if target_idx == -1:
        return []

    target_col = df_user_id.columns[target_idx]
    for value in df_user_id[target_col].astype(str).str.strip():
        if value:
            user_ids.add(value)
    return sorted(user_ids)


def read_data_sheet(excel_file):
    return pd.read_excel(
        excel_file,
        sheet_name=SHEET_NAME,
        engine="openpyxl",
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )


def display_file_label(excel_file, idx: int) -> str:
    return f"{excel_file.name} [upload {idx + 1}]"


def df_to_tsv(df: pd.DataFrame) -> str:
    if df is None or df.empty:
        return "<empty>"
    return df.to_csv(sep="\t", index=False)


def preview_tsv(df: pd.DataFrame, n: int = 20) -> str:
    if df is None or df.empty:
        return "<empty>"
    return df.head(n).to_csv(sep="\t", index=False)


def render_result_df(df: pd.DataFrame, max_rows: int = 20):
    if df is None or df.empty:
        st.caption("No rows.")
        return
    if st.session_state.get("display_as_table", False):
        st.dataframe(df.head(max_rows), width="stretch", hide_index=True)
    else:
        st.code(preview_tsv(df, n=max_rows), language="text")


def render_tsv_preview(title: str, df: pd.DataFrame, max_rows: int = 20):
    st.subheader(title)
    render_result_df(df, max_rows=max_rows)


def render_tsv_download(label: str, df: pd.DataFrame, file_name: str):
    if df is None or df.empty:
        return
    st.download_button(
        label,
        data=df_to_tsv(df).encode("utf-8"),
        file_name=file_name,
        mime="text/tab-separated-values",
    )


def is_effectively_blank_for_activity(value) -> bool:
    s = clean_cell_minimal(value)
    return s in {"", "0"}


def row_is_effectively_blank(row: pd.Series, columns) -> bool:
    for col in columns:
        if col in row.index and not is_effectively_blank_for_activity(row.get(col, "")):
            return False
    return True


def add_issue(issue_rows: list, file_name: str, row_num: str, sample_id: str, severity: str, category: str, message: str):
    issue_rows.append({
        "File": file_name,
        "Row": row_num,
        "Sample_ID": sample_id,
        "Severity": severity,
        "Category": category,
        "Message": message,
    })


def analyze_duplicates(active_df: pd.DataFrame, file_name: str, issue_rows: list):
    details = []
    if active_df.empty or "index" not in active_df.columns or "index2" not in active_df.columns:
        return details, set()

    working = active_df.copy()
    working["_index_key"] = working["index"].apply(clean_cell_minimal)
    working["_index2_key"] = working["index2"].apply(clean_cell_minimal)
    valid_mask = (
        working["_index_key"].map(lambda v: not is_effectively_blank_for_activity(v))
        & working["_index2_key"].map(lambda v: not is_effectively_blank_for_activity(v))
    )
    dup_rows = working.loc[valid_mask].loc[
        working.loc[valid_mask].duplicated(subset=["_index_key", "_index2_key"], keep=False)
    ]

    combos = set(tuple(x) for x in working.loc[valid_mask, ["_index_key", "_index2_key"]].drop_duplicates().itertuples(index=False, name=None))

    if dup_rows.empty:
        return details, combos

    for (idx1, idx2), grp in dup_rows.groupby(["_index_key", "_index2_key"], sort=False):
        rows = grp["_row_num"].astype(str).tolist()
        samples = [s for s in grp["Sample_ID"].astype(str).tolist() if s]
        message = f"Duplicate index/index2 combination {idx1}/{idx2} found in rows: {', '.join(rows)}"
        for _, rec in grp.iterrows():
            add_issue(
                issue_rows,
                file_name=file_name,
                row_num=str(rec["_row_num"]),
                sample_id=clean_cell_minimal(rec.get("Sample_ID", "")),
                severity="Error",
                category="Index Duplication",
                message=message,
            )
        details.append({
            "File": file_name,
            "index": idx1,
            "index2": idx2,
            "Rows": ", ".join(rows),
            "Sample_IDs": ", ".join(samples),
        })
    return details, combos


def cross_file_duplicates(file_combo_map: dict):
    combo_files_map = {}
    for file_label, combos in file_combo_map.items():
        for combo in combos:
            combo_files_map.setdefault(combo, []).append(file_label)
    return {combo: files for combo, files in combo_files_map.items() if len(files) > 1}


def qc_screen_file(excel_file, upload_idx: int):
    file_name = display_file_label(excel_file, upload_idx)
    user_ids = read_user_ids(excel_file)
    try:
        df_raw = read_data_sheet(excel_file)
    except ValueError:
        return {
            "file_name": file_name,
            "sheet_found": False,
            "user_ids": user_ids,
            "issues": [{
                "File": file_name,
                "Row": "",
                "Sample_ID": "",
                "Severity": "Error",
                "Category": "Input / Structural QC",
                "Message": f'Missing required sheet "{SHEET_NAME}".',
            }],
            "amplicon_summary": [],
            "duplicate_details": [],
            "row_count_checked": 0,
            "combos": set(),
            "missing_required_columns": CORE_COLS.copy(),
        }

    df = normalize_header_columns(blankify_minimal(df_raw))
    if BASE_EDITING_COL not in df.columns:
        df[BASE_EDITING_COL] = ""
    if "Exon" not in df.columns:
        df["Exon"] = ""
    if "Expected_HDR_Amplicon" not in df.columns:
        df["Expected_HDR_Amplicon"] = ""
    if NGRNA_FINAL_HEADER not in df.columns:
        df[NGRNA_FINAL_HEADER] = ""
    df["_row_num"] = df.index + 2

    issues = []
    missing_required_columns = [col for col in CORE_COLS if col not in df.columns]
    for col in missing_required_columns:
        add_issue(issues, file_name, "", "", "Error", "Input / Structural QC", f'Missing required column "{col}".')

    rows_checked = 0
    unique_amplicons = []
    amp_seen = {}

    available_relevant_cols = [c for c in QC_ANCHOR_COLS if c in df.columns]
    active_df = df.loc[~df.apply(lambda row: row_is_effectively_blank(row, available_relevant_cols), axis=1)].copy()

    for _, row in active_df.iterrows():
        rows_checked += 1
        row_num = str(int(row["_row_num"]))
        sample_id = clean_cell_minimal(row.get("Sample_ID", ""))

        missing_fields = []
        for col in CORE_COLS:
            if col not in df.columns:
                continue
            if clean_cell_minimal(row.get(col, "")) == "":
                missing_fields.append(col)
        if missing_fields:
            add_issue(issues, file_name, row_num, sample_id, "Error", "Core Row QC", f"Missing required value(s): {', '.join(missing_fields)}")

        amp_raw = clean_cell_minimal(row.get("Amplicon", ""))
        amp_norm = normalize_dna(amp_raw)
        if amp_norm:
            if amp_norm not in amp_seen:
                amp_seen[amp_norm] = f"Amplicon{len(amp_seen)+1:02d}"
                unique_amplicons.append({
                    "Amplicon Label": amp_seen[amp_norm],
                    "Length": len(amp_norm),
                    "Occurrence Count": 1,
                    "Sequence": amp_norm,
                })
            else:
                for item in unique_amplicons:
                    if item["Sequence"] == amp_norm:
                        item["Occurrence Count"] += 1
                        break

            if "U" in amp_norm:
                add_issue(issues, file_name, row_num, sample_id, "Error", "Amplicon QC", "Amplicon contains U. Amplicon must be a DNA sequence with A/T/C/G only.")
            elif not DNA_ONLY_RE.fullmatch(amp_norm):
                add_issue(issues, file_name, row_num, sample_id, "Error", "Amplicon QC", "Amplicon contains characters other than A/T/C/G.")

        grna_raw = clean_cell_minimal(row.get("gRNA", ""))
        grna_norm = normalize_dna(grna_raw)
        if grna_norm:
            if "U" in grna_norm:
                add_issue(issues, file_name, row_num, sample_id, "Info", "gRNA QC", "gRNA contains U. This will be converted to T in File Merge mode.")
            elif not DNA_ONLY_RE.fullmatch(grna_norm) and not DNA_OR_U_RE.fullmatch(grna_norm):
                add_issue(issues, file_name, row_num, sample_id, "Error", "gRNA QC", "gRNA contains characters other than A/T/C/G/U.")

            if amp_norm and DNA_ONLY_RE.fullmatch(amp_norm) and DNA_OR_U_RE.fullmatch(grna_norm):
                grna_for_match = grna_norm.replace("U", "T")
                grna_rc = reverse_complement(grna_for_match)
                if grna_for_match not in amp_norm and grna_rc not in amp_norm:
                    add_issue(issues, file_name, row_num, sample_id, "Error", "Amplicon QC", "gRNA was not found in Amplicon, and reverse complement(gRNA) was also not found.")

        exon_raw = clean_cell_minimal(row.get("Exon", ""))
        exon_norm = normalize_dna(exon_raw)
        if exon_norm:
            if "U" in exon_norm:
                add_issue(issues, file_name, row_num, sample_id, "Error", "Exon QC", "Exon contains U. Exon must be a DNA sequence with A/T/C/G only.")
            elif not DNA_ONLY_RE.fullmatch(exon_norm):
                add_issue(issues, file_name, row_num, sample_id, "Error", "Exon QC", "Exon contains characters other than A/T/C/G.")
            elif amp_norm and DNA_ONLY_RE.fullmatch(amp_norm) and exon_norm not in amp_norm:
                add_issue(issues, file_name, row_num, sample_id, "Error", "Exon QC", "Exon was not found in Amplicon in the same orientation.")

        hdr_amp_raw = clean_cell_minimal(row.get("Expected_HDR_Amplicon", ""))
        hdr_amp_norm = normalize_dna(hdr_amp_raw)
        if hdr_amp_norm:
            if "U" in hdr_amp_norm:
                add_issue(issues, file_name, row_num, sample_id, "Error", "Logical Consistency QC", "Expected_HDR_Amplicon contains U. Expected_HDR_Amplicon must be a DNA sequence with A/T/C/G only.")
            elif not DNA_ONLY_RE.fullmatch(hdr_amp_norm):
                add_issue(issues, file_name, row_num, sample_id, "Error", "Logical Consistency QC", "Expected_HDR_Amplicon contains characters other than A/T/C/G.")
            if clean_cell_minimal(row.get("Quantification_Window_Coordinates", "")) == "":
                add_issue(issues, file_name, row_num, sample_id, "Error", "Logical Consistency QC", "Expected_HDR_Amplicon is filled but Quantification_Window_Coordinates is blank.")

        base_edit_val = normalize_base_editing_type(row.get(BASE_EDITING_COL, ""))
        if base_edit_val:
            if base_edit_val not in VALID_BASE_EDITING_TYPES:
                add_issue(issues, file_name, row_num, sample_id, "Error", "Base Editing QC", f'Invalid Base_Editing_Type "{base_edit_val}". Allowed values are ABE, CBE, BOTH.')
            else:
                has_conflict = any(clean_cell_minimal(row.get(col, "")) != "" for col in BASE_EDITING_CONFLICT_COLS if col in df.columns)
                if has_conflict:
                    add_issue(issues, file_name, row_num, sample_id, "Error", "Base Editing QC", "Base-editing row also contains HDR/window/ngRNA fields.")

    dup_details, combos = analyze_duplicates(active_df, file_name, issues)

    return {
        "file_name": file_name,
        "sheet_found": True,
        "user_ids": user_ids,
        "issues": issues,
        "amplicon_summary": unique_amplicons,
        "duplicate_details": dup_details,
        "row_count_checked": rows_checked,
        "combos": combos,
        "missing_required_columns": missing_required_columns,
    }


def run_qc(uploaded_files):
    all_issues = []
    amplicon_rows = []
    duplicate_rows = []
    user_ids = set()
    file_summaries = []
    file_combo_map = {}
    rows_checked_total = 0

    for upload_idx, excel_file in enumerate(uploaded_files):
        result = qc_screen_file(excel_file, upload_idx)
        file_name = result["file_name"]
        rows_checked_total += result["row_count_checked"]
        user_ids.update(result["user_ids"])
        file_combo_map[file_name] = result["combos"]
        all_issues.extend(result["issues"])
        for row in result["amplicon_summary"]:
            amplicon_rows.append({"File": file_name, **row})
        duplicate_rows.extend(result["duplicate_details"])

        issue_df = pd.DataFrame(result["issues"])
        errors = int((issue_df["Severity"] == "Error").sum()) if not issue_df.empty else 0
        infos = int((issue_df["Severity"] == "Info").sum()) if not issue_df.empty else 0
        warnings = int((issue_df["Severity"] == "Warning").sum()) if not issue_df.empty else 0
        file_summaries.append({
            "File": file_name,
            "Rows Checked": result["row_count_checked"],
            "Errors": errors,
            "Warnings": warnings,
            "Info": infos,
            "Status": "Needs Correction" if errors > 0 else "Pass",
        })

    cross_dup = cross_file_duplicates(file_combo_map)
    for (idx1, idx2), files in cross_dup.items():
        message = f"Cross-file duplicate index/index2 combination {idx1}/{idx2} found in files: {', '.join(files)}"
        all_issues.append({
            "File": ", ".join(files),
            "Row": "",
            "Sample_ID": "",
            "Severity": "Error",
            "Category": "Index Duplication",
            "Message": message,
        })
        duplicate_rows.append({
            "File": ", ".join(files),
            "index": idx1,
            "index2": idx2,
            "Rows": "",
            "Sample_IDs": "",
        })

    issues_df = pd.DataFrame(all_issues)
    if issues_df.empty:
        issues_df = pd.DataFrame(columns=["File", "Row", "Sample_ID", "Severity", "Category", "Message"])
    amplicon_df = pd.DataFrame(amplicon_rows)
    duplicate_df = pd.DataFrame(duplicate_rows)
    summary_df = pd.DataFrame(file_summaries)

    return {
        "issues_df": issues_df,
        "amplicon_df": amplicon_df,
        "duplicate_df": duplicate_df,
        "summary_df": summary_df,
        "user_ids": sorted(user_ids),
        "files_checked": len(uploaded_files),
        "rows_checked": rows_checked_total,
        "error_count": int((issues_df["Severity"] == "Error").sum()),
        "warning_count": int((issues_df["Severity"] == "Warning").sum()),
        "info_count": int((issues_df["Severity"] == "Info").sum()),
        "cross_dup_combos": cross_dup,
    }


def clean_sample_id(value) -> str:
    s = clean_cell_merge(value)
    if not s:
        return ""
    return re.sub(r"[^A-Za-z0-9]", "", s)


def make_unique_sample_id(base_id: str, occurrence_num: int) -> str:
    if occurrence_num <= 1 or base_id == "":
        return base_id
    return f"{base_id}_{occurrence_num}"


def apply_sample_id_cleanup_and_duplicates(df: pd.DataFrame, auto_fix_duplicates: bool):
    df_out = df.copy()
    if "Sample_ID" not in df_out.columns:
        df_out["Sample_ID"] = ""

    original_ids = df_out["Sample_ID"].apply(clean_cell_merge)
    cleaned_ids = original_ids.apply(clean_sample_id)
    df_out["Sample_ID"] = cleaned_ids

    collisions_mask = cleaned_ids.ne("") & cleaned_ids.duplicated(keep=False)
    preview_rows = []
    collision_count = int(collisions_mask.sum())

    final_ids = cleaned_ids.copy()
    if auto_fix_duplicates:
        occurrence_tracker = {}
        for idx, cleaned_id in cleaned_ids.items():
            if cleaned_id == "":
                continue
            occurrence_tracker[cleaned_id] = occurrence_tracker.get(cleaned_id, 0) + 1
            final_ids.at[idx] = make_unique_sample_id(cleaned_id, occurrence_tracker[cleaned_id])
        df_out["Sample_ID"] = final_ids

    if collision_count > 0:
        for idx in df_out.index[collisions_mask].tolist():
            preview_rows.append({
                "File": clean_cell_merge(df_out.at[idx, "_source_file"]) if "_source_file" in df_out.columns else "",
                "Row": int(df_out.at[idx, "_source_row"]) if "_source_row" in df_out.columns else "",
                "Original Sample_ID": original_ids.at[idx],
                "Cleaned Sample_ID": cleaned_ids.at[idx],
                "Final Sample_ID": final_ids.at[idx],
            })

    cleanup_changed_count = int((original_ids != cleaned_ids).sum())

    return df_out, {
        "cleanup_changed_count": cleanup_changed_count,
        "collision_row_count": collision_count,
        "collision_preview_df": pd.DataFrame(preview_rows),
        "auto_fix_applied": auto_fix_duplicates and collision_count > 0,
    }


def clean_and_qc_grna_merge(df: pd.DataFrame):
    if "gRNA" not in df.columns:
        return df, {"u_to_t": 0, "invalid_after": 0, "invalid_examples": [], "missing_grna": len(df), "missing_grna_examples": []}
    grna_orig = df["gRNA"].copy()
    contains_u_mask = grna_orig.astype(str).str.contains(r"[Uu]", regex=True, na=False)
    missing_mask = grna_orig.astype(str).str.strip() == ""
    df["gRNA"] = grna_orig.apply(normalize_grna_for_merge)
    nonempty_mask = df["gRNA"].astype(str).str.len() > 0
    invalid_mask = nonempty_mask & ~df["gRNA"].astype(str).str.match(DNA_ONLY_RE)

    missing_examples = []
    if missing_mask.any():
        example_df = df.loc[missing_mask, ["Sample_ID"]].copy() if "Sample_ID" in df.columns else pd.DataFrame(index=df.index[missing_mask])
        if "_source_row" in df.columns:
            example_df["_source_row"] = df.loc[missing_mask, "_source_row"]
        for idx, row in example_df.head(5).iterrows():
            sample = clean_cell_merge(row.get("Sample_ID", ""))
            source_row = row.get("_source_row", None)
            if sample:
                missing_examples.append(sample)
            elif pd.notna(source_row):
                missing_examples.append(f"row {int(source_row)}")
            else:
                missing_examples.append(f"row {idx + 2}")

    return df, {
        "u_to_t": int(contains_u_mask.sum()),
        "invalid_after": int(invalid_mask.sum()),
        "invalid_examples": df.loc[invalid_mask, "gRNA"].astype(str).unique().tolist()[:5],
        "missing_grna": int(missing_mask.sum()),
        "missing_grna_examples": missing_examples,
    }


def apply_base_editing_rules_merge(df: pd.DataFrame):
    df_out = df.copy()
    if BASE_EDITING_COL not in df_out.columns:
        df_out[BASE_EDITING_COL] = ""
    if BE_Q30_CUTOFF_COL not in df_out.columns:
        df_out[BE_Q30_CUTOFF_COL] = ""
    df_out[BASE_EDITING_COL] = df_out[BASE_EDITING_COL].apply(lambda x: normalize_base_editing_type(clean_cell_merge(x)))

    amp_invalid_examples = []
    hdr_invalid_examples = []
    grna_not_found_examples = []
    base_edit_invalid_examples = []
    base_edit_conflict_examples = []
    be_q30_invalid_examples = []
    reoriented_examples = []

    amp_invalid_count = 0
    hdr_invalid_count = 0
    grna_not_found_count = 0
    base_edit_invalid_count = 0
    base_edit_conflict_count = 0
    be_q30_invalid_count = 0
    reoriented_count = 0

    for idx, row in df_out.iterrows():
        sample_id = clean_cell_merge(row.get("Sample_ID", "")) or f"row {int(row.get('_source_row', idx + 2))}"
        amp_raw = clean_cell_merge(row.get("Amplicon", ""))
        amp_norm = amp_raw.upper().replace("U", "T") if amp_raw else ""
        grna_raw = clean_cell_merge(row.get("gRNA", ""))
        grna_norm = normalize_grna_for_merge(grna_raw) if grna_raw else ""
        hdr_amp_raw = clean_cell_merge(row.get("Expected_HDR_Amplicon", ""))
        hdr_amp_norm = hdr_amp_raw.upper() if hdr_amp_raw else ""
        bet = normalize_base_editing_type(row.get(BASE_EDITING_COL, ""))
        be_q30_raw = clean_cell_merge(row.get(BE_Q30_CUTOFF_COL, ""))
        be_q30_norm = normalize_be_q30_cutoff(be_q30_raw)

        if amp_norm and not DNA_ONLY_RE.fullmatch(amp_norm):
            amp_invalid_count += 1
            if len(amp_invalid_examples) < 5:
                amp_invalid_examples.append(f"{sample_id}: {amp_raw}")

        if hdr_amp_norm and not DNA_ONLY_RE.fullmatch(hdr_amp_norm):
            hdr_invalid_count += 1
            if len(hdr_invalid_examples) < 5:
                hdr_invalid_examples.append(f"{sample_id}: {hdr_amp_raw}")

        if grna_norm:
            grna_rc = reverse_complement(grna_norm)
            grna_in_amp = bool(amp_norm) and (grna_norm in amp_norm)
            grna_rc_in_amp = bool(amp_norm) and (grna_rc in amp_norm)
            if amp_norm and not (grna_in_amp or grna_rc_in_amp):
                grna_not_found_count += 1
                if len(grna_not_found_examples) < 5:
                    grna_not_found_examples.append(sample_id)
            if bet in VALID_BASE_EDITING_TYPES and amp_norm:
                amp_rc = reverse_complement(amp_norm)
                grna_in_amp_rc = grna_norm in amp_rc
                if (not grna_in_amp) and grna_in_amp_rc:
                    df_out.at[idx, "Amplicon"] = amp_rc
                    reoriented_count += 1
                    if len(reoriented_examples) < 5:
                        reoriented_examples.append(sample_id)

        if bet:
            if bet not in VALID_BASE_EDITING_TYPES:
                base_edit_invalid_count += 1
                if len(base_edit_invalid_examples) < 5:
                    base_edit_invalid_examples.append(f"{sample_id}: {bet}")
            else:
                has_conflict = any(clean_cell_merge(row.get(col, "")) != "" for col in BASE_EDITING_CONFLICT_COLS if col in df_out.columns)
                if has_conflict:
                    base_edit_conflict_count += 1
                    if len(base_edit_conflict_examples) < 5:
                        base_edit_conflict_examples.append(sample_id)

        if be_q30_norm:
            if be_q30_norm not in VALID_BE_Q30_CUTOFF_VALUES:
                be_q30_invalid_count += 1
                df_out.at[idx, BE_Q30_CUTOFF_COL] = ""
                if len(be_q30_invalid_examples) < 5:
                    be_q30_invalid_examples.append(f"{sample_id}: {be_q30_raw}")
            else:
                df_out.at[idx, BE_Q30_CUTOFF_COL] = be_q30_norm
        else:
            df_out.at[idx, BE_Q30_CUTOFF_COL] = ""

    return blankify_merge(df_out), {
        "amp_invalid": amp_invalid_count,
        "amp_invalid_examples": amp_invalid_examples,
        "hdr_invalid": hdr_invalid_count,
        "hdr_invalid_examples": hdr_invalid_examples,
        "grna_not_found": grna_not_found_count,
        "grna_not_found_examples": grna_not_found_examples,
        "base_edit_invalid": base_edit_invalid_count,
        "base_edit_invalid_examples": base_edit_invalid_examples,
        "base_edit_conflict": base_edit_conflict_count,
        "base_edit_conflict_examples": base_edit_conflict_examples,
        "be_q30_invalid": be_q30_invalid_count,
        "be_q30_invalid_examples": be_q30_invalid_examples,
        "amplicon_reoriented": reoriented_count,
        "amplicon_reoriented_examples": reoriented_examples,
    }



def parse_quantification_coordinates(value):
    """Parse an inclusive 0-based coordinate string such as ``123-230``."""
    text = clean_cell_merge(value).replace("–", "-")
    match = re.fullmatch(r"(\d+)\s*-\s*(\d+)", text)
    if not match:
        return None
    return int(match.group(1)), int(match.group(2))


def find_guide_hit_for_merge(reference: str, guide: str):
    """Return guide orientation and reference positions for merge-time calculations."""
    reference_u = clean_cell_merge(reference).upper().replace("U", "T")
    guide_u = normalize_grna_for_merge(guide)
    if not reference_u or not guide_u:
        return None

    start = reference_u.find(guide_u)
    if start >= 0:
        return {
            "strand": "+",
            "start": start,
            "end": start + len(guide_u) - 1,
            "ref_5p": start,
            "ref_3p": start + len(guide_u) - 1,
        }

    guide_rc = reverse_complement(guide_u)
    start = reference_u.find(guide_rc)
    if start >= 0:
        return {
            "strand": "-",
            "start": start,
            "end": start + len(guide_u) - 1,
            "ref_5p": start + len(guide_u) - 1,
            "ref_3p": start,
        }
    return None


def blank_window_required_read_length_merge(peg_hit: dict, ref_len: int):
    """Return the minimum safe read length for a blank-window row.

    Blank quantification-window rows are treated as base-editing/cutting rows
    for this preflight only. Their implicit CRISPResso center is -10 and their
    implicit plotting half-window is 22. The signed center is converted to the
    reference axis using the gRNA strand, then the requested exclusion and
    manual safety margin are added exactly as specified by the user.
    """
    if not peg_hit:
        return None

    quant_center = BLANK_WINDOW_QUANTIFICATION_CENTER
    if peg_hit["strand"] == "+":
        absolute_center = int(peg_hit["ref_3p"] + quant_center)
    else:
        # A negative guide-relative offset moves toward increasing reference
        # coordinates on the reverse strand.
        absolute_center = int(peg_hit["ref_3p"] - quant_center)

    required_read_length = int(
        absolute_center
        + BLANK_WINDOW_PLOT_SIZE
        + BLANK_WINDOW_CRISPRESSO_EXCLUSION
        + PLOT_WINDOW_SAFETY_MARGIN
    )
    return {
        "required_read_length": required_read_length,
        "absolute_center": absolute_center,
        "plot_window_size": BLANK_WINDOW_PLOT_SIZE,
        "quantification_window_center": quant_center,
        "guide_strand": peg_hit["strand"],
        "guide_3p": int(peg_hit["ref_3p"]),
        "reference_length": int(ref_len),
    }


def signed_center_from_absolute_merge(center_position: int, peg_hit: dict, ng_hit: dict, fallback_center: int) -> int:
    """Convert an absolute plotting center to CRISPResso's signed guide-relative offset."""
    if ng_hit is not None and ng_hit["ref_5p"] != peg_hit["ref_5p"]:
        direction = 1 if ng_hit["ref_5p"] > peg_hit["ref_5p"] else -1
    elif fallback_center != 0:
        strand_flip = 1 if peg_hit["strand"] == "+" else -1
        direction = 1 if fallback_center * strand_flip > 0 else -1
    else:
        direction = 0

    if direction == 0:
        return 0
    strand_flip = 1 if peg_hit["strand"] == "+" else -1
    distance = abs(int(center_position) - peg_hit["ref_3p"])
    return int(distance * direction * strand_flip)


def recenter_plot_parameters_merge(
    plot_window_size: int,
    quant_center: int,
    peg_hit: dict,
    ng_hit: dict,
    ref_len: int,
):
    """Keep or recenter a symmetric plot interval against the full [0, ref_len) reference."""
    half = max(0, int(plot_window_size))
    direction = 0
    if ng_hit is not None and ng_hit["ref_5p"] != peg_hit["ref_5p"]:
        direction = 1 if ng_hit["ref_5p"] > peg_hit["ref_5p"] else -1
    elif quant_center != 0:
        strand_flip = 1 if peg_hit["strand"] == "+" else -1
        direction = 1 if quant_center * strand_flip > 0 else -1

    original_center = int(peg_hit["ref_3p"] + direction * abs(int(quant_center)))
    left = original_center - half
    right = original_center + half

    # Use strict interior bounds.  In particular, CRISPResso can reject an
    # exact right-boundary case (e.g. 160 + 55 == ref_len == 215).
    if left > 0 and right < ref_len:
        return int(quant_center), int(plot_window_size), False, ""

    clipped_left = max(0, left)
    clipped_right = min(ref_len, right)
    if clipped_left >= clipped_right:
        new_center = max(0, min(original_center, max(0, ref_len - 1)))
        new_half = 0
    else:
        new_center = (clipped_left + clipped_right) // 2
        boundary_limited_half = min(
            new_center - clipped_left,
            clipped_right - new_center,
        )
        # Reduce the deduced boundary-limited half-window by 2 bp for safety.
        new_half = max(0, boundary_limited_half - PLOT_WINDOW_SAFETY_MARGIN)

    new_quant_center = signed_center_from_absolute_merge(
        new_center, peg_hit, ng_hit, quant_center
    )
    message = (
        f"Plot window [{left}, {right}) reached or exceeded full trimmed "
        f"Amplicon bounds [0, {ref_len}); recentered at {new_center} and "
        f"applied a {PLOT_WINDOW_SAFETY_MARGIN}-bp safety margin, changing "
        f"Quantification_Window_Center to {new_quant_center} and "
        f"Plot_Window_Size to {new_half}."
    )
    return int(new_quant_center), int(new_half), True, message


def apply_read_length_rules_merge(df: pd.DataFrame, read_length: int):
    """Apply read-length trimming and make CRISPResso windows safe for each row."""
    df_out = df.copy()
    read_length = max(0, int(read_length))
    stats = {
        "amplicon_trimmed": 0,
        "amplicon_bases_removed": 0,
        "hdr_trimmed": 0,
        "exon_adjusted": 0,
        "quant_coordinates_clamped": 0,
        "plot_recentered": 0,
        "plot_recenter_examples": [],
        "blank_window_rows_checked": 0,
        "blank_window_read_length_danger": 0,
        "blank_window_danger_examples": [],
        "grna_not_found_after_trim": 0,
        "grna_not_found_after_trim_examples": [],
    }

    for idx, row in df_out.iterrows():
        amp_raw = clean_cell_merge(row.get("Amplicon", ""))
        hdr_raw = clean_cell_merge(row.get("Expected_HDR_Amplicon", ""))
        amp_original = amp_raw
        original_amp_len = len(amp_original)
        trim_count = max(0, original_amp_len - read_length)
        amp_trimmed = amp_original[:read_length] if trim_count else amp_original

        if trim_count:
            df_out.at[idx, "Amplicon"] = amp_trimmed
            stats["amplicon_trimmed"] += 1
            stats["amplicon_bases_removed"] += trim_count
            if hdr_raw:
                # Hard rule: Expected_HDR_Amplicon loses exactly the same number
                # of 3' bases as Amplicon, independent of HDR edit type.
                df_out.at[idx, "Expected_HDR_Amplicon"] = hdr_raw[:max(0, len(hdr_raw) - trim_count)]
                stats["hdr_trimmed"] += 1

            exon = clean_cell_merge(row.get("Exon", ""))
            exon_start = amp_original.upper().find(exon.upper()) if exon else -1
            if exon_start >= 0 and exon_start < len(amp_trimmed):
                exon_end = exon_start + len(exon)
                retained_exon = amp_original[exon_start:min(exon_end, len(amp_trimmed))]
                if retained_exon != exon:
                    df_out.at[idx, "Exon"] = retained_exon
                    stats["exon_adjusted"] += 1
            elif exon and exon in amp_original and len(amp_trimmed) <= exon_start:
                df_out.at[idx, "Exon"] = ""
                stats["exon_adjusted"] += 1

        amp_for_match = clean_cell_merge(df_out.at[idx, "Amplicon"])
        grna = normalize_grna_for_merge(row.get("gRNA", ""))
        if grna and amp_for_match:
            grna_rc = reverse_complement(grna)
            if grna not in amp_for_match.upper() and grna_rc not in amp_for_match.upper():
                stats["grna_not_found_after_trim"] += 1
                if len(stats["grna_not_found_after_trim_examples"]) < 5:
                    stats["grna_not_found_after_trim_examples"].append(
                        clean_cell_merge(row.get("Sample_ID", "")) or f"row {int(row.get('_source_row', idx + 2))}"
                    )

        coords = parse_quantification_coordinates(row.get("Quantification_Window_Coordinates", ""))
        if coords and amp_for_match:
            start, end = coords
            ref_len = len(amp_for_match)
            new_start = max(0, min(start, ref_len - 1))
            # Quantification_Window_Coordinates are inclusive and 0-based.
            # Keep the start unchanged (apart from ordinary reference bounds),
            # but exclude the final 15 reference bases from the end coordinate.
            eligible_right = max(0, ref_len - CRISPRESSO_EXCLUDE_BP_FROM_RIGHT - 1)
            new_end = max(new_start, min(end, eligible_right))
            new_coords = f"{new_start}-{new_end}"
            if new_coords != clean_cell_merge(row.get("Quantification_Window_Coordinates", "")):
                df_out.at[idx, "Quantification_Window_Coordinates"] = new_coords
                stats["quant_coordinates_clamped"] += 1

        qcenter_text = clean_cell_merge(row.get("Quantification_Window_Center", ""))
        plot_text = clean_cell_merge(row.get("Plot_Window_Size", ""))
        peg_hit = find_guide_hit_for_merge(amp_original, grna)

        # Base-editing and cutting rows intentionally leave all three window
        # fields blank. Do not populate those fields here; instead, preflight
        # the selected read length using the agreed implicit -10 / 22 values.
        # The calculation is based on the original amplicon so a guide that is
        # removed by 3' trimming can still be evaluated correctly.
        blank_window_row = not qcenter_text and not plot_text and not clean_cell_merge(
            row.get("Quantification_Window_Coordinates", "")
        )
        if blank_window_row:
            stats["blank_window_rows_checked"] += 1
            blank_check = blank_window_required_read_length_merge(peg_hit, original_amp_len)
            if blank_check is not None:
                effective_read_length = len(amp_for_match)
                if effective_read_length < blank_check["required_read_length"]:
                    stats["blank_window_read_length_danger"] += 1
                    if len(stats["blank_window_danger_examples"]) < 5:
                        sample = clean_cell_merge(row.get("Sample_ID", "")) or f"row {int(row.get('_source_row', idx + 2))}"
                        stats["blank_window_danger_examples"].append(
                            f"{sample}: read length {read_length} (effective {effective_read_length}) < "
                            f"minimum {blank_check['required_read_length']}; "
                            f"center {blank_check['absolute_center']}, "
                            f"plot window {blank_check['plot_window_size']}, "
                            f"gRNA {blank_check['guide_strand']} strand 3' end {blank_check['guide_3p']}"
                        )
            continue

        try:
            qcenter = int(float(qcenter_text))
            plot_size = int(float(plot_text))
        except (TypeError, ValueError):
            continue

        if not peg_hit or not amp_for_match:
            continue
        ng_hit = find_guide_hit_for_merge(amp_original, row.get(NGRNA_FINAL_HEADER, ""))
        new_qcenter, new_plot_size, recentered, message = recenter_plot_parameters_merge(
            plot_size, qcenter, peg_hit, ng_hit, len(amp_for_match)
        )
        if recentered:
            df_out.at[idx, "Quantification_Window_Center"] = str(new_qcenter)
            df_out.at[idx, "Plot_Window_Size"] = str(new_plot_size)
            stats["plot_recentered"] += 1
            if len(stats["plot_recenter_examples"]) < 5:
                sample = clean_cell_merge(row.get("Sample_ID", "")) or f"row {int(row.get('_source_row', idx + 2))}"
                stats["plot_recenter_examples"].append(f"{sample}: {message}")

    return blankify_merge(df_out), stats

def collect_removed_miseq_rows(df: pd.DataFrame):
    removed_records = []
    valid_mask = []
    for _, row in df.iterrows():
        missing_cols = [col for col in MISEQ_REQUIRED_COLS if clean_cell_merge(row.get(col, "")) == ""]
        is_valid = len(missing_cols) == 0
        valid_mask.append(is_valid)
        if not is_valid:
            source_row = row.get("_source_row", "")
            source_file = clean_cell_merge(row.get("_source_file", ""))
            sample_id = clean_cell_merge(row.get("Sample_ID", ""))
            removed_records.append({
                "File": source_file,
                "Row": int(source_row) if str(source_row) != "" else "",
                "Sample_ID": sample_id,
                "Missing Fields": ", ".join(missing_cols),
            })
    return pd.Series(valid_mask, index=df.index), pd.DataFrame(removed_records)


def run_merge(uploaded_files, auto_fix_duplicates=False):
    raw_list, exp_list, logs = [], [], []
    file_combos = {}
    grna_qc_logs = []
    base_edit_qc_logs = []
    user_ids = set()
    sample_id_preview_rows = []

    total_u_to_t = 0
    total_invalid = 0
    total_missing_grna = 0
    total_invalid_examples = set()
    total_missing_grna_examples = []
    total_amp_invalid = 0
    total_hdr_invalid = 0
    total_grna_not_found = 0
    total_base_edit_invalid = 0
    total_base_edit_conflict = 0
    total_be_q30_invalid = 0
    total_amplicon_reoriented = 0
    total_amplicon_trimmed = 0
    total_amplicon_bases_removed = 0
    total_hdr_trimmed = 0
    total_exon_adjusted = 0
    total_quant_coordinates_clamped = 0
    total_plot_recentered = 0
    total_plot_recenter_examples = []
    total_blank_window_rows_checked = 0
    total_blank_window_read_length_danger = 0
    total_blank_window_danger_examples = []
    total_grna_not_found_after_trim = 0
    total_grna_not_found_after_trim_examples = []
    total_amp_invalid_examples = set()
    total_hdr_invalid_examples = set()
    total_grna_not_found_examples = set()
    total_base_edit_invalid_examples = set()
    total_base_edit_conflict_examples = set()
    total_be_q30_invalid_examples = set()
    total_amplicon_reoriented_examples = set()
    sample_id_cleanup_changed_count = 0
    sample_id_collision_row_count = 0
    sample_id_auto_fix_applied = False

    for upload_idx, excel_file in enumerate(uploaded_files):
        file_label = display_file_label(excel_file, upload_idx)
        try:
            df_data = read_data_sheet(excel_file)
        except ValueError:
            logs.append({
                "File": file_label, "Sheet Found": False, "Input Samples": 0,
                "Missing gRNA": 0, "In-file Dup Combos": 0, "Cross-file Dup": False,
                "gRNA U→T": 0, "Non-ATCG (post U→T)": 0,
                "Amplicon Non-ATCG": 0, "Expected_HDR_Amplicon Non-ATCG": 0,
                "gRNA Not Found in Amplicon": 0, "Invalid Base_Editing_Type": 0,
                "Base Editing Conflicts": 0, "Invalid BE_Q30_cutoff": 0, "Amplicons Reoriented": 0,
            })
            continue

        user_ids.update(read_user_ids(excel_file))
        df = normalize_header_columns(blankify_merge(df_data))
        if BASE_EDITING_COL not in df.columns:
            df[BASE_EDITING_COL] = ""
        if BE_Q30_CUTOFF_COL not in df.columns:
            df[BE_Q30_CUTOFF_COL] = ""
        if "Exon" not in df.columns:
            df["Exon"] = ""
        if "Expected_HDR_Amplicon" not in df.columns:
            df["Expected_HDR_Amplicon"] = ""
        if NGRNA_FINAL_HEADER not in df.columns:
            df[NGRNA_FINAL_HEADER] = ""

        missing_core_cols = [col for col in CORE_COLS if col not in df.columns]
        for col in missing_core_cols:
            df[col] = ""

        df["_source_row"] = df.index + 2
        df["_source_file"] = file_label

        relevant_cols = [c for c in CORE_COLS + QC_OPTIONAL_COLS if c in df.columns]
        mask = pd.Series(False, index=df.index)
        if relevant_cols:
            mask = df[relevant_cols].apply(lambda r: any(clean_cell_merge(v) != "" for v in r), axis=1)
        df_nonblank = blankify_merge(df[mask].copy())

        df_nonblank, qc = clean_and_qc_grna_merge(df_nonblank)
        raw_list.append(df_nonblank.copy())

        dup_details, combos = analyze_duplicates(df_nonblank.assign(_row_num=df_nonblank["_source_row"]), file_label, [])
        file_combos[file_label] = combos
        in_dup_count = len(dup_details)

        df_processed, qc_be = apply_base_editing_rules_merge(df_nonblank)
        df_processed, trim_qc = apply_read_length_rules_merge(
            df_processed, int(st.session_state.get("miseq_reads", 284))
        )
        exp_list.append(df_processed)

        total_u_to_t += qc["u_to_t"]
        total_invalid += qc["invalid_after"]
        total_missing_grna += qc["missing_grna"]
        total_invalid_examples.update(qc["invalid_examples"])
        for ex in qc["missing_grna_examples"]:
            if ex not in total_missing_grna_examples and len(total_missing_grna_examples) < 5:
                total_missing_grna_examples.append(ex)
        total_amp_invalid += qc_be["amp_invalid"]
        total_amp_invalid_examples.update(qc_be["amp_invalid_examples"])
        total_hdr_invalid += qc_be["hdr_invalid"]
        total_hdr_invalid_examples.update(qc_be["hdr_invalid_examples"])
        total_grna_not_found += qc_be["grna_not_found"]
        total_grna_not_found_examples.update(qc_be["grna_not_found_examples"])
        total_base_edit_invalid += qc_be["base_edit_invalid"]
        total_base_edit_invalid_examples.update(qc_be["base_edit_invalid_examples"])
        total_base_edit_conflict += qc_be["base_edit_conflict"]
        total_base_edit_conflict_examples.update(qc_be["base_edit_conflict_examples"])
        total_be_q30_invalid += qc_be["be_q30_invalid"]
        total_be_q30_invalid_examples.update(qc_be["be_q30_invalid_examples"])
        total_amplicon_reoriented += qc_be["amplicon_reoriented"]
        total_amplicon_reoriented_examples.update(qc_be["amplicon_reoriented_examples"])
        total_amplicon_trimmed += trim_qc["amplicon_trimmed"]
        total_amplicon_bases_removed += trim_qc["amplicon_bases_removed"]
        total_hdr_trimmed += trim_qc["hdr_trimmed"]
        total_exon_adjusted += trim_qc["exon_adjusted"]
        total_quant_coordinates_clamped += trim_qc["quant_coordinates_clamped"]
        total_plot_recentered += trim_qc["plot_recentered"]
        for example in trim_qc["plot_recenter_examples"]:
            if example not in total_plot_recenter_examples and len(total_plot_recenter_examples) < 5:
                total_plot_recenter_examples.append(example)
        total_blank_window_rows_checked += trim_qc["blank_window_rows_checked"]
        total_blank_window_read_length_danger += trim_qc["blank_window_read_length_danger"]
        for example in trim_qc["blank_window_danger_examples"]:
            if example not in total_blank_window_danger_examples and len(total_blank_window_danger_examples) < 5:
                total_blank_window_danger_examples.append(example)
        total_grna_not_found_after_trim += trim_qc["grna_not_found_after_trim"]
        for example in trim_qc["grna_not_found_after_trim_examples"]:
            if example not in total_grna_not_found_after_trim_examples and len(total_grna_not_found_after_trim_examples) < 5:
                total_grna_not_found_after_trim_examples.append(example)

        grna_qc_logs.append({
            "File": file_label,
            "Missing gRNA": qc["missing_grna"],
            "gRNA U→T": qc["u_to_t"],
            "Non-ATCG (post U→T)": qc["invalid_after"],
            "Examples (up to 5)": ", ".join(qc["invalid_examples"]) if qc["invalid_examples"] else "",
        })
        base_edit_qc_logs.append({
            "File": file_label,
            "Amplicon Non-ATCG": qc_be["amp_invalid"],
            "Expected_HDR_Amplicon Non-ATCG": qc_be["hdr_invalid"],
            "gRNA Not Found in Amplicon": qc_be["grna_not_found"],
            "Invalid Base_Editing_Type": qc_be["base_edit_invalid"],
            "Base Editing Conflicts": qc_be["base_edit_conflict"],
            "Invalid BE_Q30_cutoff": qc_be["be_q30_invalid"],
            "Amplicons Reoriented": qc_be["amplicon_reoriented"],
            "gRNA Not Found After Trimming": trim_qc["grna_not_found_after_trim"],
            "Amplicons Trimmed": trim_qc["amplicon_trimmed"],
            "HDR Amplicons Trimmed": trim_qc["hdr_trimmed"],
            "Exons Adjusted": trim_qc["exon_adjusted"],
            "Quantification Coordinates Clamped": trim_qc["quant_coordinates_clamped"],
            "Plot Windows Recentered": trim_qc["plot_recentered"],
            "Blank Window Rows Checked": trim_qc["blank_window_rows_checked"],
            "Blank Window Read-Length Danger": trim_qc["blank_window_read_length_danger"],
        })
        logs.append({
            "File": file_label, "Sheet Found": True, "Input Samples": len(df_nonblank),
            "Missing gRNA": qc["missing_grna"], "In-file Dup Combos": in_dup_count, "Cross-file Dup": False,
            "gRNA U→T": qc["u_to_t"], "Non-ATCG (post U→T)": qc["invalid_after"],
            "Amplicon Non-ATCG": qc_be["amp_invalid"],
            "Expected_HDR_Amplicon Non-ATCG": qc_be["hdr_invalid"],
            "gRNA Not Found in Amplicon": qc_be["grna_not_found"],
            "Invalid Base_Editing_Type": qc_be["base_edit_invalid"],
            "Base Editing Conflicts": qc_be["base_edit_conflict"],
            "Invalid BE_Q30_cutoff": qc_be["be_q30_invalid"],
            "Amplicons Reoriented": qc_be["amplicon_reoriented"],
            "gRNA Not Found After Trimming": trim_qc["grna_not_found_after_trim"],
            "Amplicons Trimmed": trim_qc["amplicon_trimmed"],
            "HDR Amplicons Trimmed": trim_qc["hdr_trimmed"],
            "Exons Adjusted": trim_qc["exon_adjusted"],
            "Quantification Coordinates Clamped": trim_qc["quant_coordinates_clamped"],
            "Plot Windows Recentered": trim_qc["plot_recentered"],
            "Blank Window Rows Checked": trim_qc["blank_window_rows_checked"],
            "Blank Window Read-Length Danger": trim_qc["blank_window_read_length_danger"],
        })

    cross_dup = cross_file_duplicates(file_combos)
    for entry in logs:
        entry["Cross-file Dup"] = any(entry["File"] in files for files in cross_dup.values())

    raw_df = blankify_merge(pd.concat(raw_list, ignore_index=True)) if raw_list else pd.DataFrame(columns=FINAL_COLS)
    expanded_df = blankify_merge(pd.concat(exp_list, ignore_index=True)) if exp_list else pd.DataFrame(columns=FINAL_COLS)

    # Sample_ID cleanup + duplicate resolution pooled across ALL uploaded files (global scope).
    expanded_df, sample_id_meta = apply_sample_id_cleanup_and_duplicates(expanded_df, auto_fix_duplicates)
    sample_id_cleanup_changed_count = sample_id_meta["cleanup_changed_count"]
    sample_id_collision_row_count = sample_id_meta["collision_row_count"]
    sample_id_auto_fix_applied = sample_id_meta["auto_fix_applied"]
    if not raw_df.empty and "Sample_ID" in expanded_df.columns and len(raw_df) == len(expanded_df):
        raw_df["Sample_ID"] = expanded_df["Sample_ID"].values

    valid_mask, removed_miseq_df = collect_removed_miseq_rows(expanded_df)
    output_df = blankify_merge(expanded_df.loc[valid_mask].copy()) if not expanded_df.empty else expanded_df.copy()
    raw_valid_df = blankify_merge(raw_df.loc[valid_mask].copy()) if not raw_df.empty else raw_df.copy()

    log_df = pd.DataFrame(logs)
    if not log_df.empty:
        totals = {"File": "Total"}
        for col in log_df.columns:
            if col not in ("File", "Sheet Found", "Cross-file Dup"):
                totals[col] = int(pd.to_numeric(log_df[col], errors="coerce").fillna(0).sum())
        log_rows = pd.concat([log_df, pd.DataFrame([totals])], ignore_index=True)
    else:
        log_rows = pd.DataFrame(columns=["File"])

    collision_preview_df = sample_id_meta["collision_preview_df"]

    return {
        "raw_df": raw_valid_df,
        "expanded_df": output_df,
        "log_rows": log_rows,
        "grna_qc_logs": pd.DataFrame(grna_qc_logs) if grna_qc_logs else pd.DataFrame(),
        "base_edit_qc_logs": pd.DataFrame(base_edit_qc_logs) if base_edit_qc_logs else pd.DataFrame(),
        "grna_qc_totals": {
            "missing_grna": total_missing_grna,
            "missing_grna_examples": total_missing_grna_examples[:5],
            "u_to_t": total_u_to_t,
            "invalid_after": total_invalid,
            "invalid_examples": sorted(list(total_invalid_examples))[:5],
        },
        "base_edit_qc_totals": {
            "amp_invalid": total_amp_invalid,
            "amp_invalid_examples": sorted(list(total_amp_invalid_examples))[:5],
            "hdr_invalid": total_hdr_invalid,
            "hdr_invalid_examples": sorted(list(total_hdr_invalid_examples))[:5],
            "grna_not_found": total_grna_not_found,
            "grna_not_found_examples": sorted(list(total_grna_not_found_examples))[:5],
            "base_edit_invalid": total_base_edit_invalid,
            "base_edit_invalid_examples": sorted(list(total_base_edit_invalid_examples))[:5],
            "base_edit_conflict": total_base_edit_conflict,
            "base_edit_conflict_examples": sorted(list(total_base_edit_conflict_examples))[:5],
            "be_q30_invalid": total_be_q30_invalid,
            "be_q30_invalid_examples": sorted(list(total_be_q30_invalid_examples))[:5],
            "amplicon_reoriented": total_amplicon_reoriented,
            "amplicon_reoriented_examples": sorted(list(total_amplicon_reoriented_examples))[:5],
            "amplicon_trimmed": total_amplicon_trimmed,
            "amplicon_bases_removed": total_amplicon_bases_removed,
            "hdr_trimmed": total_hdr_trimmed,
            "exon_adjusted": total_exon_adjusted,
            "quant_coordinates_clamped": total_quant_coordinates_clamped,
            "plot_recentered": total_plot_recentered,
            "plot_recenter_examples": total_plot_recenter_examples[:5],
            "blank_window_rows_checked": total_blank_window_rows_checked,
            "blank_window_read_length_danger": total_blank_window_read_length_danger,
            "blank_window_danger_examples": total_blank_window_danger_examples[:5],
            "grna_not_found_after_trim": total_grna_not_found_after_trim,
            "grna_not_found_after_trim_examples": total_grna_not_found_after_trim_examples[:5],
        },
        "cross_dup_combos": cross_dup,
        "user_ids": sorted(user_ids),
        "sample_id_cleanup_changed_count": sample_id_cleanup_changed_count,
        "sample_id_collision_row_count": sample_id_collision_row_count,
        "sample_id_collision_preview_df": collision_preview_df,
        "sample_id_auto_fix_applied": sample_id_auto_fix_applied,
        "removed_miseq_df": removed_miseq_df,
    }


def render_qc_results(results):
    st.subheader("QC Summary")
    st.write(
        f"Files Checked: {results['files_checked']} | "
        f"Rows Checked: {results['rows_checked']} | "
        f"Errors: {results['error_count']} | "
        f"Warnings: {results['warning_count']} | "
        f"Info: {results['info_count']}"
    )

    if results["user_ids"]:
        st.info(f"Unique User IDs found: {', '.join(results['user_ids'])}")

    if not results["summary_df"].empty:
        render_tsv_preview("Per-file Status", results["summary_df"], max_rows=200)
        render_tsv_download("Download per-file status TSV", results["summary_df"], "qc_per_file_status.tsv")

    if not results["issues_df"].empty:
        st.subheader("QC Issues")
        show_non_errors = st.checkbox("Show Warnings and Info", value=False)
        issues_df = results["issues_df"].copy()
        if not show_non_errors:
            issues_df = issues_df[issues_df["Severity"] == "Error"].copy()
        issues_df["_severity_sort"] = issues_df["Severity"].map(SEVERITY_ORDER)
        issues_df = issues_df.sort_values(by=["_severity_sort", "File", "Row", "Category", "Message"], kind="stable")
        issues_df["Severity"] = issues_df["Severity"].map(SEVERITY_LABELS).fillna(issues_df["Severity"])
        issues_df = issues_df.drop(columns=["_severity_sort"])
        render_result_df(issues_df, max_rows=200)
        render_tsv_download("Download QC issues TSV", issues_df, "qc_issues.tsv")
    else:
        st.success("No QC issues found.")

    if not results["amplicon_df"].empty:
        display_cols = ["File", "Amplicon Label", "Length", "Occurrence Count", "Sequence"]
        amp_df = results["amplicon_df"][display_cols].copy()
        render_tsv_preview("Unique Amplicon Length Summary", amp_df, max_rows=200)
        render_tsv_download("Download amplicon summary TSV", amp_df, "qc_amplicon_summary.tsv")

    if not results["duplicate_df"].empty:
        render_tsv_preview("Duplicate Index/index2 Details", results["duplicate_df"], max_rows=200)
        render_tsv_download("Download duplicate details TSV", results["duplicate_df"], "qc_duplicate_details.tsv")

def render_merge_results(results, prefix):
    final_df = ensure_columns_and_order(results["expanded_df"], FINAL_COLS)

    if results["sample_id_collision_row_count"] > 0:
        if results["sample_id_auto_fix_applied"]:
            st.warning(
                f"⚠️Found **{results['sample_id_collision_row_count']}** row(s) with duplicate cleaned Sample_ID values. "
                f"Auto-correction was applied because the duplicate-cleanup option is enabled."
            )
        else:
            st.warning(
                f"⚠️Found **{results['sample_id_collision_row_count']}** row(s) with duplicate cleaned Sample_ID values after removing non-alphanumeric characters. "
                f"Enable the duplicate-cleanup option if you want the tool to auto-correct them."
            )
        if not results["sample_id_collision_preview_df"].empty:
            render_tsv_preview("Sample_ID Collision Preview", results["sample_id_collision_preview_df"], max_rows=200)
            render_tsv_download("Download Sample_ID collision TSV", results["sample_id_collision_preview_df"], "sample_id_collisions.tsv")

    grna_totals = results["grna_qc_totals"]
    if grna_totals["missing_grna"] > 0:
        examples = ", ".join(grna_totals["missing_grna_examples"])
        msg = f"⚠️Found **{grna_totals['missing_grna']}** row(s) with missing gRNA. gRNA is required for CRISPResso."
        if examples:
            msg += f" Examples: {examples}"
        st.warning(msg)
    if grna_totals["u_to_t"] > 0:
        st.info(f"gRNA cleanup: Converted **{grna_totals['u_to_t']}** entries from U→T (case-insensitive).")
    if grna_totals["invalid_after"] > 0:
        examples = ", ".join(grna_totals["invalid_examples"])
        msg = f"Found **{grna_totals['invalid_after']}** gRNA entry(ies) with non-A/T/C/G characters after U→T conversion."
        if examples:
            msg += f" Examples: {examples}"
        st.warning(msg)

    be_totals = results["base_edit_qc_totals"]
    if be_totals["amp_invalid"] > 0:
        examples = ", ".join(be_totals["amp_invalid_examples"])
        msg = f"⚠️Found **{be_totals['amp_invalid']}** Amplicon entry(ies) with non-A/T/C/G characters."
        if examples:
            msg += f" Examples: {examples}"
        st.warning(msg)
    if be_totals["hdr_invalid"] > 0:
        examples = ", ".join(be_totals["hdr_invalid_examples"])
        msg = f"⚠️Found **{be_totals['hdr_invalid']}** Expected_HDR_Amplicon entry(ies) with non-A/T/C/G characters."
        if examples:
            msg += f" Examples: {examples}"
        st.warning(msg)
    if be_totals["grna_not_found"] > 0:
        examples = ", ".join(be_totals["grna_not_found_examples"])
        msg = f"⚠️Found **{be_totals['grna_not_found']}** row(s) where gRNA was not found in Amplicon and reverse complement(gRNA) was also not found."
        if examples:
            msg += f" Examples: {examples}"
        st.warning(msg)
    if be_totals["base_edit_invalid"] > 0:
        examples = ", ".join(be_totals["base_edit_invalid_examples"])
        msg = f"⚠️Found **{be_totals['base_edit_invalid']}** row(s) with invalid Base_Editing_Type. Allowed values are ABE, CBE, or BOTH."
        if examples:
            msg += f" Examples: {examples}"
        st.warning(msg)
    if be_totals["base_edit_conflict"] > 0:
        examples = ", ".join(be_totals["base_edit_conflict_examples"])
        msg = f"⚠️Found **{be_totals['base_edit_conflict']}** base-editing row(s) that also contain HDR/window/ngRNA fields."
        if examples:
            msg += f" Examples: {examples}"
        st.warning(msg)
    if be_totals["be_q30_invalid"] > 0:
        examples = ", ".join(be_totals["be_q30_invalid_examples"])
        msg = f"⚠️Found **{be_totals['be_q30_invalid']}** row(s) with invalid BE_Q30_cutoff values. Allowed values are ON or OFF; invalid entries were cleared to blank."
        if examples:
            msg += f" Examples: {examples}"
        st.warning(msg)
    if be_totals["amplicon_reoriented"] > 0:
        examples = ", ".join(be_totals["amplicon_reoriented_examples"])
        msg = f"Re-oriented **{be_totals['amplicon_reoriented']}** Amplicon row(s) so gRNA is in the same sense for base-editing samples."
        if examples:
            msg += f" Examples: {examples}"
        st.info(msg)
    if be_totals["amplicon_trimmed"] > 0:
        st.info(
            f"Read-length trimming: trimmed **{be_totals['amplicon_trimmed']}** Amplicon row(s), "
            f"removing **{be_totals['amplicon_bases_removed']}** total 3′ bases; "
            f"Expected_HDR_Amplicon was trimmed by the same per-row amount."
        )
    if be_totals["exon_adjusted"] > 0:
        st.info(f"Adjusted **{be_totals['exon_adjusted']}** Exon value(s) to remain within the retained Amplicon sequence.")
    if be_totals["quant_coordinates_clamped"] > 0:
        st.info(f"Clamped **{be_totals['quant_coordinates_clamped']}** quantification coordinate range(s) to the retained Amplicon.")
    if be_totals["plot_recentered"] > 0:
        examples = "; ".join(be_totals["plot_recenter_examples"])
        msg = f"Recentered **{be_totals['plot_recentered']}** plot window(s) against the full trimmed Amplicon [0, L)."
        if examples:
            msg += f" Examples: {examples}"
        st.info(msg)
    if be_totals["blank_window_read_length_danger"] > 0:
        examples = "; ".join(be_totals["blank_window_danger_examples"])
        msg = (
            f"🚨 DANGER: **{be_totals['blank_window_read_length_danger']}** base-editing/cutting row(s) "
            f"have blank quantification-window fields, but the selected read length is too short for "
            f"the implicit center (-10), plot window (22), 5-bp CRISPResso exclusion, and "
            f"{PLOT_WINDOW_SAFETY_MARGIN}-bp safety margin. The resulting CRISPResso plot may fail or omit "
            f"part of the intended window."
        )
        if examples:
            msg += f" Examples: {examples}"
        st.error(msg)

    if be_totals["grna_not_found_after_trim"] > 0:
        examples = ", ".join(be_totals["grna_not_found_after_trim_examples"])
        msg = f"⚠️After trimming, the full gRNA was not found in **{be_totals['grna_not_found_after_trim']}** Amplicon row(s)."
        if examples:
            msg += f" Examples: {examples}"
        st.warning(msg)

    if not results["log_rows"].empty and "In-file Dup Combos" in results["log_rows"].columns:
        per_file_dup_df = results["log_rows"][results["log_rows"]["File"] != "Total"].copy()
        dup_mask = pd.to_numeric(per_file_dup_df["In-file Dup Combos"], errors="coerce").fillna(0) > 0
        if dup_mask.any():
            dup_rows = per_file_dup_df.loc[dup_mask, ["File", "In-file Dup Combos"]]
            details = ", ".join(
                f"{row['File']} ({int(pd.to_numeric(row['In-file Dup Combos'], errors='coerce'))})"
                for _, row in dup_rows.head(5).iterrows()
            )
            total_dup_files = int(dup_mask.sum())
            st.warning(f"⚠️Found in-file duplicate index/index2 combos in **{total_dup_files}** file(s). Examples: {details}")

    if results["cross_dup_combos"]:
        combos_str = ", ".join(f"{i}/{j}" for i, j in results["cross_dup_combos"])
        st.warning(f"⚠️Cross-file duplicate index combos: {combos_str}")

    hdr_mask = final_df["Expected_HDR_Amplicon"].astype(str).str.strip().ne("") if not final_df.empty else pd.Series(dtype=bool)
    qwc_blank = final_df["Quantification_Window_Coordinates"].astype(str).str.strip().eq("") if not final_df.empty else pd.Series(dtype=bool)
    prob_mask = hdr_mask & qwc_blank
    if prob_mask.any():
        n_prob = int(prob_mask.sum())
        examples = final_df.loc[prob_mask, "Sample_ID"].astype(str).head(5).tolist()
        st.warning(f"⚠️{n_prob} row(s) have Expected_HDR_Amplicon but missing Quantification_Window_Coordinates. Examples: {', '.join(examples)}")

    if not results["removed_miseq_df"].empty:
        n_removed = len(results["removed_miseq_df"])
        preview = results["removed_miseq_df"].head(5).copy()
        examples = "; ".join(
            f"{row['File']} row {row['Row']} missing {row['Missing Fields']}"
            for _, row in preview.iterrows()
        )
        st.warning(
            f"⚠️Removed **{n_removed}** row(s) from the merge output because MiSeq requires Sample_ID, I7_Index_ID, index, I5_Index_ID, and index2. "
            f"Examples: {examples}"
        )
        with st.expander("Show removed rows for MiSeq-required fields"):
            render_result_df(results["removed_miseq_df"], max_rows=200)
            render_tsv_download("Download removed MiSeq rows TSV", results["removed_miseq_df"], "removed_miseq_rows.tsv")

    if results["user_ids"]:
        st.info(f"Merged Unique User IDs: {', '.join(results['user_ids'])}")

    render_tsv_preview("Merge Log", results["log_rows"], max_rows=200)
    render_tsv_download("Download merge log TSV", results["log_rows"], "merge_log.tsv")

    if not results["grna_qc_logs"].empty:
        with st.expander("Show per-file gRNA QC details"):
            render_result_df(results["grna_qc_logs"], max_rows=200)
            render_tsv_download("Download gRNA QC TSV", results["grna_qc_logs"], "grna_qc.tsv")
    if not results["base_edit_qc_logs"].empty:
        with st.expander("Show per-file Base Editing / Amplicon QC details"):
            render_result_df(results["base_edit_qc_logs"], max_rows=200)
            render_tsv_download("Download base editing QC TSV", results["base_edit_qc_logs"], "base_edit_qc.tsv")

    excel_buf = BytesIO()
    write_excel_with_blanks(final_df, excel_buf, SHEET_NAME, results["user_ids"])
    st.download_button(
        "📥 Download merged Excel (CRISPResso)",
        data=excel_buf,
        file_name=out_excel,
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

    csv_buf = StringIO()
    csv_buf.write("[Header]\n")
    csv_buf.write(f"Experiment Name,{prefix}\n")
    csv_buf.write(f"Date,{today.month}/{today.day}/{today.year}\n")
    csv_buf.write("Workflow,GenerateFASTQ\n")
    csv_buf.write(f"[Reads]\n{int(st.session_state.get('miseq_reads', 284))}\n")
    csv_buf.write("[Settings]\n")
    csv_buf.write("[Data]\n")
    csv_buf.write("Sample_ID,Sample_Name,I7_Index_ID,index,I5_INDEX_ID,index2,Sample_Project,Description\n")
    for _, row in results["raw_df"].iterrows():
        csv_buf.write(
            f"{row['Sample_ID']},,{row['I7_Index_ID']},{row['index']},{row['I5_Index_ID']},{row['index2']},,\n"
        )
    st.download_button(
        "📥 Download Miseq CSV",
        data=csv_buf.getvalue().encode("utf-8"),
        file_name=out_csv,
        mime="text/csv",
    )

    render_tsv_preview("Merged Data Preview", final_df, max_rows=200)
    render_tsv_download("Download merged data TSV", final_df, f"{prefix}_merged_data.tsv")



if run_clicked:
    if not uploaded_files:
        st.warning("⚠️Please upload at least one file before running.")
    elif mode == "QC Screening":
        state.qc_results = run_qc(uploaded_files)
        state.merge_results = None
    else:
        state.merge_results = run_merge(uploaded_files, auto_fix_duplicates=auto_fix_sample_id_dups)
        state.qc_results = None

if mode == "QC Screening" and state.qc_results is not None:
    render_qc_results(state.qc_results)
elif mode == "File Merge" and state.merge_results is not None:
    render_merge_results(state.merge_results, prefix)
