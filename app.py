import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

st.set_page_config(layout="wide")
st.title("🧬 DFU Biofilm: Empirical Inhibition Engine")

@st.cache_data
def load_data():
    return pd.read_csv("all_gene_counts_2.tsv", sep='\t')

df = load_data()

# 1. EMPIRICAL INPUTS
with st.sidebar:
    st.header("Laboratory Parameters")
    # Control vs Treated counts
    ctrl_count = st.number_input("Baseline Bacterial Count (CFU/mL)", value=1e6)
    treat_count = st.number_input("Treated Bacterial Count (CFU/mL)", value=1e4)
    gene = st.selectbox("Target Gene", df['gene_name'].unique())

# 2. THE REAL BIOFILM CALCULATION (Log-Reduction)
if st.button("Calculate Log-Reduction"):
    # The standard microbiological formula: Log10(Control) - Log10(Treated)
    log_reduction = np.log10(ctrl_count) - np.log10(treat_count)
    
    st.subheader(f"Log-Reduction Value: {log_reduction:.2f}")
    st.write("Interpretation: A value > 3 is typically considered bactericidal for biofilm inhibition.")

    # 3. EMPIRICAL VISUALIZATION
    # This plots your actual data (Control vs Treatment)
    st.subheader("Biofilm Suppression Analysis")
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.barplot(x=['Baseline', 'Treated'], y=[ctrl_count, treat_count], palette="Blues")
    ax.set_yscale('log')
    ax.set_ylabel("CFU/mL (Log Scale)")
    st.pyplot(fig)
