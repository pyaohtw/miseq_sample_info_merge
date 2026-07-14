# MiSeq Sample Info Merger

A Streamlit-based utility for **QC screening** and **MiSeq / CRISPResso file preparation** from one or more Excel workbooks.

The app supports two main workflows:

- **QC Screening**: read-only validation of uploaded sample sheets
- **File Merge**: merge uploaded sheets into CRISPResso-ready Excel output and MiSeq-ready CSV output

## What the App Does

### QC Screening

QC Screening checks uploaded `.xlsx` files and reports issues without modifying the source files.

It currently checks for:

- Missing required `Data` sheet
- Missing required columns
- Missing required values in core columns
- Duplicate `index/index2` combinations within a file
- Duplicate `index/index2` combinations across files
- Invalid or inconsistent DNA sequence fields
- `gRNA` not found in `Amplicon` (including reverse-complement check)
- `Exon` not found in `Amplicon`
- Invalid `Base_Editing_Type`
- Base-editing rows that also contain conflicting HDR/window/ngRNA fields
- `Expected_HDR_Amplicon` filled while `Quantification_Window_Coordinates` is blank

QC output is shown in a **safe text/TSV-style preview** by default for better stability in hosted environments.

### File Merge

File Merge combines uploaded `Data` sheets into two outputs:

1. **Merged Excel (CRISPResso)**
   - Sheet name: `Data`
   - Columns normalized into the app's final output schema
   - Blank-like values cleaned up
   - Optional `UserID` sheet preserved as a merged unique list

2. **MiSeq CSV**
   - Includes the standard MiSeq header sections
   - Uses only rows with all required MiSeq fields present
   - Uses the configured single-end read length

Additional merge behavior includes:

- Normalizing `ngRNA` header variants into `ngRNA`
- Cleaning `Sample_ID` by removing non-alphanumeric characters
- Optional auto-correction of duplicate cleaned `Sample_ID` values
- Converting `gRNA` values from `U` to `T`
- Base-editing QC checks during merge
- Re-orienting `Amplicon` for base-editing rows when needed
- Removing rows from MiSeq CSV output if required MiSeq fields are missing

## Current Features

1. **Two operating modes**
   - `QC Screening`
   - `File Merge`

2. **Multiple workbook upload**
   Upload one or more `.xlsx` files in a single run.

3. **Safe result rendering**
   Result previews are displayed as TSV-style text blocks instead of interactive dataframes/tables by default.

4. **Downloadable outputs**
   - Merged CRISPResso Excel file
   - MiSeq CSV file
   - TSV downloads for QC summaries, issue lists, and merge logs

5. **MiSeq read length control**
   - Default single-end read length: **284**
   - Adjustable in File Merge mode

6. **Date-based default naming**
   Default file prefix format:

   ```text
   runYYMMDDmi
   ```

## Input Expectations

### Required workbook structure

Each uploaded workbook should contain:

- A sheet named `Data`
- Optionally, a sheet named `UserID`

### Core required columns

The app expects these core columns for full QC / merge behavior:

```text
Sample_ID
I7_Index_ID
index
I5_Index_ID
index2
Amplicon
gRNA
```

### MiSeq-required columns

Rows must contain these values to be included in the MiSeq CSV:

```text
Sample_ID
I7_Index_ID
index
I5_Index_ID
index2
```

## Getting Started

### Prerequisites

- Python 3.10+
- Conda or another virtual environment manager
- Git

### Installation

1. Clone the repo:

   ```bash
   git clone https://github.com/your-org/miseq-sample-info-merger.git
   cd miseq-sample-info-merger
   ```

2. Create and activate an environment:

   ```bash
   conda create -n miseq python=3.10
   conda activate miseq
   ```

3. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

### Recommended requirements

If you are running into hosted-environment instability, use pinned dependencies in `requirements.txt`.

Example stable set:

```txt
streamlit==1.39.0
pandas==2.2.2
numpy==1.26.4
openpyxl==3.1.5
xlrd==2.0.1
pytz==2024.1
```

### Optional Streamlit config

For some hosted environments, the following file can improve stability:

```text
.streamlit/config.toml
```

```toml
[server]
fileWatcherType = "none"
```

## Usage

Run the app with:

```bash
streamlit run streamlit_app.py
```

Then open the local URL shown in the console, usually:

```text
http://localhost:8501
```

## Repo Structure

```text
.
├── streamlit_app.py    # Main Streamlit app
├── requirements.txt    # Python dependencies
├── .streamlit/
│   └── config.toml     # Optional Streamlit runtime config
└── README.md           # This file
```

## Output Summary

### QC Screening outputs

- QC summary
- Per-file status
- QC issues list
- Unique amplicon summary
- Duplicate index/index2 details
- TSV-style previews and TSV downloads for result tables

### File Merge outputs

- Merged Excel workbook for CRISPResso
- MiSeq sample sheet CSV
- Merge log
- gRNA QC summary
- Base editing / amplicon QC summary
- TSV-style preview of merged data

## Notes

- The app treats blank-like values such as `nan`, `none`, `null`, `0`, and `0.0` as empty in several cleaning steps.
- `gRNA` values containing `U` are converted to `T` during merge.
- The app preserves text formatting in the merged Excel output where possible.
- Hosted environments may behave differently from local or Codespaces runs, especially around Streamlit frontend rendering.
