import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats

st.title("🧬 DFU Biofilm Research Engine: Validated Statistical Pipeline")

@st.cache_data
def load_and_normalize():
    df = pd.read_csv("all_gene_counts_2.tsv", sep='\t')
    
    # Identify sample columns
    sample_cols = [c for c in df.columns if c not in ['gene_id', 'gene_name']]
    
    # NORMALIZATION: CPM (Counts Per Million)
    # This adjusts for varying sequencing depth between X1, X2, etc.
    df_norm = df.copy()
    for col in sample_cols:
        df_norm[col] = (df[col] / df[col].sum()) * 1e6
        
    return df_norm, sample_cols

df, sample_cols = load_and_normalize()

# 2. DEFINING GROUPS (Replace with your specific design)
control_cols = ['X1', 'X2', 'X3'] 
treatment_cols = ['X4', 'X5', 'X6']

# 3. STATISTICAL ENGINE
# Log2 Fold Change on normalized data
df['log2FC'] = np.log2((df[treatment_cols].mean(axis=1) + 1) / (df[control_cols].mean(axis=1) + 1))

# T-test
def calculate_p(row):
    _, p = stats.ttest_ind(row[control_cols], row[treatment_cols], equal_var=False)
    return p

df['p_value'] = df.apply(calculate_p, axis=1)

# -log10 transformation for Volcano Plot
df['sig_metric'] = -np.log10(df['p_value'])

# 4. DASHBOARD
st.write("### Statistically Validated Results")
st.dataframe(df[['gene_name', 'log2FC', 'sig_metric']].sort_values('sig_metric', ascending=False))

# Visualization
st.write("### Volcano Plot (Industry Standard)")
st.scatter_chart(df, x='log2FC', y='sig_metric')
