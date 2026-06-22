import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats

st.title("🧬 DFU Biofilm Research Engine: Robust Analysis")

@st.cache_data
def load_and_analyze():
    df = pd.read_csv("all_gene_counts_2.tsv", sep='\t')
    
    # 1. Clean data: Ensure columns are numbers
    control_cols = ['X1', 'X2', 'X3'] 
    treatment_cols = ['X4', 'X5', 'X6']
    
    for col in control_cols + treatment_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)

    # 2. Normalize: CPM
    for col in control_cols + treatment_cols:
        df[col] = (df[col] / df[col].sum()) * 1e6

    # 3. Robust Statistical Engine
    def calculate_p(row):
        try:
            # We filter out rows that are all zeros to prevent math errors
            c_data = row[control_cols]
            t_data = row[treatment_cols]
            if c_data.sum() == 0 and t_data.sum() == 0:
                return 1.0
            _, p = stats.ttest_ind(c_data, t_data, equal_var=False)
            return p
        except:
            return 1.0

    df['log2FC'] = np.log2((df[treatment_cols].mean(axis=1) + 1) / (df[control_cols].mean(axis=1) + 1))
    df['p_value'] = df.apply(calculate_p, axis=1)
    df['sig_metric'] = -np.log10(df['p_value'].replace(0, 1e-300)) # Protect against log(0)
    
    return df

df = load_and_analyze()

# 4. Dashboard
st.write("### Statistically Validated Results")
st.dataframe(df[['gene_name', 'log2FC', 'p_value']].sort_values('p_value'))

st.write("### Volcano Plot")
st.scatter_chart(df, x='log2FC', y='sig_metric')
