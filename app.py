import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from scipy import stats

st.set_page_config(layout="wide")
st.title("🧬 DFU Biofilm Research Engine: PhD Analytical Pipeline")

@st.cache_data
def load_and_analyze():
    # Load your actual research data
    df = pd.read_csv("all_gene_counts_2.tsv", sep='\t')
    return df

df = load_and_analyze()

# Sidebar Setup
st.sidebar.header("Experimental Parameters")
cols = [c for c in df.columns if c not in ['gene_id', 'gene_name']]
control_cols = st.sidebar.multiselect("Select Control Replicates", cols, default=cols[0:3])
treatment_cols = st.sidebar.multiselect("Select Treatment Replicates", cols, default=cols[3:6])

def calculate_p(row):
    try:
        # Robust safety check for T-test
        c_data = pd.to_numeric(row[control_cols], errors='coerce').fillna(0)
        t_data = pd.to_numeric(row[treatment_cols], errors='coerce').fillna(0)
        if c_data.std() == 0 and t_data.std() == 0:
            return 1.0
        _, p = stats.ttest_ind(c_data, t_data, equal_var=False)
        return p
    except:
        return 1.0

if st.sidebar.button("Run Authentic Analysis"):
    # CPM Normalization
    df_norm = df.copy()
    for col in control_cols + treatment_cols:
        df_norm[col] = (df[col] / df[col].sum()) * 1e6
        
    # Statistical Calculation
    df['log2FC'] = np.log2((df_norm[treatment_cols].mean(axis=1) + 1) / (df_norm[control_cols].mean(axis=1) + 1))
    df['p_value'] = df.apply(calculate_p, axis=1)
    df['neg_log10_p'] = -np.log10(df['p_value'].replace(0, 1e-300))
    
    # Professional Volcano Plot
    st.subheader("Volcano Plot: Interaction Assessment")
    fig = px.scatter(df, x='log2FC', y='neg_log10_p', hover_name='gene_name', 
                     color='neg_log10_p', title="Differential Gene Expression Interaction")
    st.plotly_chart(fig, use_container_width=True)
    
    st.subheader("Statistically Significant Targets")
    st.dataframe(df[df['p_value'] < 0.05].sort_values('p_value'))
