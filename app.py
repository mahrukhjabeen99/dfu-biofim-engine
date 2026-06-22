import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from scipy import stats

st.set_page_config(layout="wide")
st.title("🧬 DFU Biofilm: Authentic Diagnostic Portal")

@st.cache_data
def load_data():
    return pd.read_csv("all_gene_counts_2.tsv", sep='\t')

df = load_data()

# SIDEBAR: Clinical Metadata (The "Experimental Design")
with st.sidebar:
    st.header("Experimental Mapping")
    poly = st.selectbox("Polymer Type", ["Chitosan", "Zeta", "PGL"])
    wound = st.selectbox("Wound Type", ["DFU", "Pressure"])
    
    # User-defined columns to compare
    ctrl_cols = st.multiselect("Select Control Replicates", df.columns[2:], default=[df.columns[2]])
    target_cols = st.multiselect("Select Treated Replicates", df.columns[2:], default=[df.columns[3]])

# DIAGNOSTIC ENGINE
if st.button("Run Authentic Analysis"):
    # 1. AUDIT GATE: Check for Empty/Zero Data
    if df[ctrl_cols + target_cols].isnull().values.any():
        st.error("Audit Failed: Your data contains missing values (NaN).")
    elif df[ctrl_cols + target_cols].sum().sum() == 0:
        st.error("Audit Failed: The selected columns contain only zeros.")
    else:
        # 2. RAW DATA DISTRIBUTION (The "Truth" Check)
        st.subheader("Data Variance Audit")
        st.write("Before calculating, we visualize the raw distribution. If these boxes do not overlap, your samples are fundamentally different.")
        raw_viz = df[ctrl_cols + target_cols].melt()
        fig_raw = px.box(raw_viz, x='variable', y='value', log_y=True, title="Raw Replicate Distribution")
        st.plotly_chart(fig_raw)
        
        # 3. EMPIRICAL CALCULATION
        df['log2FC'] = np.log2((df[target_cols].mean(axis=1) + 1) / (df[ctrl_cols].mean(axis=1) + 1))
        df['p_val'] = df.apply(lambda row: stats.ttest_ind(row[ctrl_cols], row[target_cols], equal_var=False)[1], axis=1)
        
        # 4. VOLCANO PLOT (The PhD Standard)
        st.subheader("Differential Expression Landscape")
        fig_volc = px.scatter(df, x='log2FC', y=-np.log10(df['p_val'] + 1e-10), 
                             hover_name='gene_name', color=df['p_val'] < 0.05,
                             title="Statistical Significance of Gene Response")
        st.plotly_chart(fig_volc)
        st.dataframe(df[df['p_val'] < 0.05].sort_values('p_val'))
