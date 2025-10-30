'''QPCR ANALYSIS (ΔCt / ΔΔCt METHOD)'''
'''owner: Nathalie Zgoda'''
'''purpose: Template script for GP2025'''
'''this version works for 96-well qPCR plates with multiple genes'''

# ---------------------------------------------------------------------
# This script demonstrates how to analyze qPCR data in Python.
# It calculates ΔCt, ΔΔCt, and relative expression (2^-ΔΔCt)
# for multiple genes and experimental conditions.
#
# The script:
#   - Reads an Excel file with Cq (Ct) values
#   - Groups data by gene and condition
#   - Calculates ΔCt normalized to a reference gene
#   - Calculates ΔΔCt relative to a calibrator condition
#   - Outputs an Excel file with all calculated results
#   - Does NOT create any plots (table output only)
#
# Note:
#   This version is designed for reproducible qPCR analysis
#   and works well with 96-well datasets.
# ---------------------------------------------------------------------


# --- 1. IMPORT REQUIRED LIBRARIES ------------------------------------
import pandas as pd       # for reading/writing Excel and handling tables
import numpy as np        # for numerical operations and error propagation


# --- 2. USER SETTINGS -------------------------------------------------
# Adjust these values to match your experiment setup.

INPUT_FILE = 'qpcr_filled.xlsx'   # <-- your filled Excel file (based on qpcr_template.xlsx)
REF_GENE   = 'ACTB'               # <-- reference gene (e.g., ACTB, GAPDH)
CALIBRATOR = 'Control'            # <-- calibrator condition (e.g., untreated cells)
USE_CQMEAN_IF_AVAILABLE = True    # use “Cq Mean” if available, otherwise “Cq”


# --- 3. LOAD YOUR DATA -----------------------------------------------
# The Excel file should include:
# Well | Cq | Cq Mean | Gene | Condition | Replicate                        #####IMPORTANT!!!
df = pd.read_excel(INPUT_FILE)


# Clean column names (remove trailing spaces etc.)
df.columns = [str(c).strip() for c in df.columns]


# --- 4. CHOOSE WHICH Ct COLUMN TO USE --------------------------------
if USE_CQMEAN_IF_AVAILABLE and 'Cq Mean' in df.columns:
    cq_col = 'Cq Mean'
elif 'Cq' in df.columns:
    cq_col = 'Cq'
else:
    raise ValueError("No 'Cq' or 'Cq Mean' column found in the input file!")


# --- 5. CONVERT TEXT NUMBERS (COMMAS → DOTS) --------------------------
# Some qPCR exports use commas as decimal separators.
def _to_num(x):
    if isinstance(x, str):
        x = x.replace(',', '.')
    return pd.to_numeric(x, errors='coerce')

# Create a new numeric column
df['Cq_num'] = df[cq_col].apply(_to_num)


# --- 6. VALIDATE REQUIRED COLUMNS ------------------------------------
for col in ['Gene', 'Condition']:
    if col not in df.columns:
        raise ValueError(f"Missing column: {col}. Please complete your Excel file.")

# Add a Replicate column if not present (optional)
if 'Replicate' not in df.columns:
    df['Replicate'] = np.nan

# Remove rows without valid Cq values
df = df.dropna(subset=['Cq_num']).copy()


# --- 7. SUMMARIZE BY CONDITION AND GENE -------------------------------
# For each Condition × Gene pair:
#   Compute mean Ct, standard deviation, and number of replicates.
grp = df.groupby(['Condition', 'Gene'], dropna=False)
summary = grp['Cq_num'].agg(['mean', 'std', 'count']).reset_index()
summary = summary.rename(columns={'mean': 'Ct_mean', 'std': 'Ct_sd', 'count': 'n'})


# --- 8. MERGE REFERENCE GENE VALUES -----------------------------------
# For each condition, merge the mean Ct of the reference gene.
ref = summary[summary['Gene'] == REF_GENE][['Condition', 'Ct_mean', 'Ct_sd', 'n']]
ref = ref.rename(columns={'Ct_mean': 'Ref_mean', 'Ct_sd': 'Ref_sd', 'n': 'Ref_n'})

# Merge with main summary
merged = pd.merge(summary, ref, on='Condition', how='left')

# Check that every condition includes the reference gene
if merged['Ref_mean'].isna().any():
    missing = merged.loc[merged['Ref_mean'].isna(), 'Condition'].unique()
    raise ValueError(f"Missing reference gene ({REF_GENE}) for these conditions: {missing}")


# --- 9. CALCULATE ΔCt AND ITS ERROR ----------------------------------
# ΔCt = Ct_gene − Ct_ref
merged['dCt'] = merged['Ct_mean'] - merged['Ref_mean']

# Propagate the standard error:
# SEM_ΔCt = sqrt( (sd_gene² / n_gene) + (sd_ref² / n_ref) )
merged['SEM_dCt'] = np.sqrt(
    (merged['Ct_sd']**2 / merged['n']) + (merged['Ref_sd']**2 / merged['Ref_n'])
)


# --- 10. CALCULATE ΔΔCt RELATIVE TO CALIBRATOR ------------------------
# Find ΔCt values for the calibrator condition
cal = merged[merged['Condition'] == CALIBRATOR][['Gene', 'dCt', 'SEM_dCt']]
cal = cal.rename(columns={'dCt': 'dCt_cal', 'SEM_dCt': 'SEM_dCt_cal'})

# Merge back into main dataframe
merged = pd.merge(merged, cal, on='Gene', how='left')

# Check that calibrator data exist for all genes
if merged['dCt_cal'].isna().any():
    missing = merged.loc[merged['dCt_cal'].isna(), 'Gene'].unique()
    raise ValueError(f"Missing calibrator ({CALIBRATOR}) for these genes: {missing}")

# ΔΔCt = ΔCt_condition − ΔCt_calibrator
merged['ddCt'] = merged['dCt'] - merged['dCt_cal']

# Propagate SEM for ΔΔCt:
# SEM_ΔΔCt = sqrt( SEM_ΔCt_condition² + SEM_ΔCt_calibrator² )
merged['SEM_ddCt'] = np.sqrt(merged['SEM_dCt']**2 + merged['SEM_dCt_cal']**2)


# --- 11. CALCULATE RELATIVE EXPRESSION (2^-ΔΔCt) ----------------------
merged['FoldChange'] = 2.0 ** (-merged['ddCt'])

# Propagate error into linear scale (approximation):
# SE_Fold = ln(2) × 2^-ΔΔCt × SEM_ΔΔCt
merged['Fold_SE'] = np.log(2) * merged['FoldChange'] * merged['SEM_ddCt']


# --- 12. PREPARE OUTPUT TABLES ---------------------------------------
tidy_cols = [
    'Condition', 'Gene', 'Ct_mean', 'Ct_sd', 'n',
    'Ref_mean', 'Ref_sd', 'Ref_n',
    'dCt', 'SEM_dCt', 'ddCt', 'SEM_ddCt',
    'FoldChange', 'Fold_SE'
]
tidy = merged[tidy_cols].sort_values(['Gene', 'Condition']).reset_index(drop=True)

# Create a wide version (pivot: genes as rows, conditions as columns)
wide = tidy.pivot_table(index='Gene', columns='Condition', values='FoldChange', aggfunc='first')


# --- 13. WRITE RESULTS TO EXCEL --------------------------------------
with pd.ExcelWriter('qpcr_results.xlsx', engine='xlsxwriter') as writer:
    tidy.to_excel(writer, index=False, sheet_name='Results_tidy')
    wide.to_excel(writer, sheet_name='FoldChange_wide')

    legend = pd.DataFrame({
        'Field': [
            'Ct_mean', 'Ct_sd', 'n', 'Ref_mean', 'Ref_sd', 'Ref_n',
            'dCt', 'ddCt', 'SEM_dCt', 'SEM_ddCt', 'FoldChange', 'Fold_SE'
        ],
        'Description': [
            'Mean Ct per Gene×Condition',
            'Standard deviation of Ct',
            'Number of replicates',
            'Mean Ct of reference gene (same condition)',
            'SD of reference gene',
            'n of reference gene',
            'ΔCt = Ct_gene - Ct_ref',
            'ΔΔCt = ΔCt_condition - ΔCt_calibrator',
            'Error of ΔCt (propagated)',
            'Error of ΔΔCt (propagated)',
            'Relative expression (2^-ΔΔCt)',
            'Linearized error for fold change'
        ]
    })
    legend.to_excel(writer, index=False, sheet_name='Legend')

print("Analysis complete.")
