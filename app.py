import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

st.set_page_config(layout="wide")
st.title("🧬 DFU Biofilm: Clinical Research Integration Portal")

# 1. LOAD DATA & METADATA
@st.cache_data
def load_data():
    return pd.read_csv("all_gene_counts_2.tsv", sep='\t')

df = load_data()

# 2. 8-FIELD RESEARCHER INTERFACE
with st.sidebar:
    st.header("Researcher Input Fields")
    wound = st.selectbox("1. Wound Type", ["Diabetic Foot Ulcer", "Pressure Ulcer"])
    bacteria = st.selectbox("2. Bacteria", ["S. aureus", "P. aeruginosa", "Both"])
    b_count = st.number_input("3. Measured Bacterial Count", value=0.0)
    poly_type = st.selectbox("4. Polymer/Nano Type", ["Chitosan", "Zeta", "PGL"])
    poly_val = st.number_input("5. Polymer Value (Zeta Potential/Conc)", value=0.0)
    gene = st.selectbox("6. Select Target Gene", df['gene_name'].unique())
    gene_val = st.number_input("7. Gene Quantitative Value", value=0.0)

# 3. Z-SCORE ANALYSIS (Objective Statistical Result)
if st.button("Generate Authentic Assessment"):
    # Filter data for the selected gene
    target_gene_data = df[df['gene_name'] == gene].iloc[:, 2:].values.flatten()
    
    # Calculate Mean and Std Dev from REAL experimental data
    mu = np.mean(target_gene_data)
    sigma = np.std(target_gene_data)
    
    # Calculate Z-score
    z_score = (gene_val - mu) / sigma if sigma != 0 else 0
    
    st.subheader(f"Final Quantitative Assessment (Z-Score): {z_score:.4f}")
    st.write("Interpretation: Z-score > 2 indicates a statistically significant deviation from baseline.")

    # 4. PHD-LEVEL VISUALIZATION (Blue Heatmap of Interactions)
    st.subheader("Biofilm Interaction Landscape")
    
    # Integration matrix of experimental inputs
    interaction_df = pd.DataFrame({
        'Parameters': ['Bacterial Count', 'Polymer Value', 'Gene Z-Score'],
        'Observed': [b_count, poly_val, z_score]
    }).set_index('Parameters')
    
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.heatmap(interaction_df, annot=True, cmap="Blues", fmt=".2f", cbar=True)
    st.pyplot(fig)
