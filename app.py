import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

st.title("🧬 DFU Biofilm Research Engine: PhD Analytical Tool")

# 1. LOAD DATA
@st.cache_data
def load_data():
    return pd.read_csv("all_gene_counts_2.tsv", sep='\t')

df = load_data()

# 2. RESEARCHER INTERFACE (Inputs)
st.sidebar.header("Experimental Parameters")
wound_type = st.sidebar.selectbox("Select Wound Type", ["Pressure Ulcer", "Diabetic Foot Ulcer"])
bacteria = st.sidebar.selectbox("Select Bacterial Species", ["P. aeruginosa", "S. aureus", "Multispecies"])
polymer = st.sidebar.text_input("Enter Polymer Nanostructure ID", "e.g., Poly-N1")

# 3. QUANTITATIVE ANALYSIS
st.write(f"### Results for {wound_type} | {bacteria} | {polymer}")

# Logic: Filter your dataframe based on the user's research inputs
# (You will need to ensure your TSV has columns for 'Wound_Type' and 'Bacteria')
filtered_df = df # Add .query() here once you map your metadata columns

st.dataframe(filtered_df.head(10))

# 4. VISUALIZATION
st.write("### Impact on Biofilm Expression")
fig, ax = plt.subplots()
# Plotting the expression counts
sns.barplot(data=filtered_df.head(10), x='gene_name', y='X1')
plt.xticks(rotation=45)
st.pyplot(fig)
