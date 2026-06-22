import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

st.title("🧬 DFU Biofilm Research Engine: PhD Analytical Tool")

# 1. LOAD YOUR REAL DATA
@st.cache_data
def load_data():
    return pd.read_csv("all_gene_counts_2.tsv", sep='\t')

df = load_data()

# 2. RESEARCHER INTERFACE (Inputs)
st.sidebar.header("Experimental Parameters")
wound_type = st.sidebar.selectbox("Select Wound Type", ["Diabetic Foot Ulcer", "Pressure Ulcer"])
bacteria_type = st.sidebar.selectbox("Select Bacteria", ["P. aeruginosa", "S. aureus"])
polymer_type = st.sidebar.selectbox("Select Polymer", ["PGL", "Chitosan"])
polymer_dose = st.sidebar.number_input("Enter Polymer Dosage (µg/mL)", value=0.0)
input_count = st.sidebar.number_input("Enter Measured Bacterial Count", value=0)

# 3. QUANTITATIVE ANALYSIS LOGIC
if st.sidebar.button("Run Authentic Analysis"):
    st.write(f"### Results for: {wound_type} | {bacteria_type} | {polymer_type} ({polymer_dose} µg/mL)")
    
    # Calculate Fold Change relative to a control baseline (example: 1000)
    baseline = 1000 
    fold_change = np.log2((input_count + 1) / (baseline + 1))
    
    # Quantitative Metric
    col1, col2 = st.columns(2)
    col1.metric("Bacterial Count", input_count)
    col2.metric("Log2 Fold Change", round(fold_change, 2))
    
    # 4. VISUALIZATION
    st.write("### Impact on Biofilm Resilience")
    fig, ax = plt.subplots(figsize=(8, 4))
    # Visual representation: Comparison of Input vs Baseline
    sns.barplot(x=['Baseline', 'Observed'], y=[baseline, input_count], palette='viridis')
    plt.ylabel("Bacterial Counts")
    st.pyplot(fig)
else:
    st.info("Enter experimental parameters in the sidebar to generate results.")
