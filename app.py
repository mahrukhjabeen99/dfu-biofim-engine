import streamlit as st
import pandas as pd
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt

st.set_page_config(layout="wide")
st.title("🧬 DFU Biofilm: PhD Research Engine (Objective Analysis)")

@st.cache_data
def load_data():
    return pd.read_csv("all_gene_counts_2.tsv", sep='\t')

df = load_data()

# 1. THE 8 INPUTS (Mapping inputs to data variables)
with st.sidebar:
    st.header("Researcher Input Fields")
    wound = st.selectbox("1. Wound Type", ["DFU", "Pressure"])
    bacteria = st.selectbox("2. Bacteria Type", ["S. aureus", "P. aeruginosa", "Both"])
    b_count = st.number_input("3. Bacterial Load (CFU/mL)", value=1e6)
    poly_type = st.selectbox("4. Polymer Type", ["Chitosan", "Zeta", "PGL"])
    poly_val = st.number_input("5. Polymer Dosage (µg/mL)", value=0.1)
    gene = st.selectbox("6. Select Target Gene", df['gene_name'].unique())
    gene_val = st.number_input("7. Measured Gene Value", value=0.0)
    run_analysis = st.button("8. Execute Objective Regression")

# 2. STATISTICAL ANALYSIS ENGINE
if run_analysis:
    st.subheader("Objective Regression: Expression ~ Clinical Inputs")
    
    # Isolate data for the selected gene
    gene_data = df[df['gene_name'] == gene].iloc[0, 2:].values
    
    # Regression model to see if inputs significantly explain variance
    X = pd.DataFrame({'B_Load': [b_count]*len(gene_data), 'Poly_Dose': [poly_val]*len(gene_data)})
    X = sm.add_constant(X)
    model = sm.OLS(gene_data, X).fit()
    
    # Display the REAL statistical truth
    st.write(model.summary())
    
    # Visualization: Heatmap of real correlation coefficients
    st.subheader("Objective Correlation Matrix")
    data_viz = pd.concat([X[['B_Load', 'Poly_Dose']], pd.Series(gene_data, name='Expression')], axis=1)
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.heatmap(data_viz.corr(), annot=True, cmap="coolwarm", ax=ax)
    st.pyplot(fig)
