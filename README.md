# MiSeq Sample Info Merger

A Streamlit-based utility to:

* **Merge** multiple Excel files’ `Data` sheets into one consolidated workbook
* **Auto-fill** blank `gRNA` entries by mapping them to every unique guide in the same Amplicon
* **Export** both an expanded Excel file (for downstream analysis) and a lightweight CSV (for MiSeq run setup)
* **QC** on your index pairs:

  * Flags in-file duplicate `index/index2` combos
  * Alerts if the same index combo appears across different files

## Features

1. **Drag-and-drop UI**
   Upload as many `.xlsx` workbooks as you like, then click **Merge**.
2. **gRNA Expansion**
   Any row missing a `gRNA` is duplicated once per unique guide found in that Amplicon, ensuring you capture all treatments and controls.
3. **Dual Output**

   * **Excel**: full expanded table, sheet named `Data`.
   * **MiSeq CSV**: fixed header (Experiment Name, Date, Workflow, Reads, Settings), then only the core indices columns for demultiplexing.
4. **Detailed Log**
   Per-file summary including:

   * Number of input samples
   * Blank-gRNA rows found
   * gRNA entries added
   * In-file and cross-file duplicate index combos
5. **Date-based default naming**
   Defaults to `runYYMMDDmi` style prefixes (e.g. `run250718mi`), but you can override.

## Getting Started

### Prerequisites

* Python 3.10+
* [Miniconda](https://docs.conda.io) or any virtual-env manager
* Git (to clone this repo)

### Installation

1. Clone the repo:

   ```bash
   git clone https://github.com/your-org/miseq-sample-info-merger.git
   cd miseq-sample-info-merger
   ```
2. Create & activate an isolated environment:

   ```bash
   conda create -n miseq python=3.10
   conda activate miseq
   ```
3. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

### Usage

```bash
streamlit run app.py
```

* Point your browser to the URL printed in the console (usually `http://localhost:8501`).
* Upload your `.xlsx` files, click **Merge**, then use the **Download** buttons to grab your Excel and CSV outputs.

## Repo Structure

```
.
├── app.py              # Main Streamlit app
├── requirements.txt    # Dependencies: streamlit, pandas, openpyxl, xlrd
└── README.md           # This file
```

## How It Works

1. **Filtering**
   Retains only rows where these columns aren’t empty:
   `Sample_ID, I7_Index_ID, index, I5_Index_ID, index2, Amplicon`
2. **gRNA Logic**

   * **Non-blank** guides pass through unchanged.
   * **Blank** guides are replaced by one row per unique `gRNA` found in the same Amplicon.
3. **QC**

   * Counts duplicate `(index,index2)` within each file.
   * Flags if any identical index pair occurs across multiple files.
4. **Export**

   * **Excel**: all fields plus filled-in `gRNA` expansions.
   * **CSV**: minimal columns for MiSeq, preserving only the original (pre-expanded) rows.


