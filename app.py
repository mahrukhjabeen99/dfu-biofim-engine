import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

st.title("🧬 DFU Biofilm Research Engine")

# 1. Load the data
@st.cache_data
def load_data():
    return pd.read_csv("all_gene_counts_2.tsv", sep='\t')

df = load_data()

# 2. Dynamic Selection based on your ACTUAL columns
st.sidebar.header("Filter Results")
# We use the actual column names from your file for selection
sample_cols = [c for c in df.columns if c not in ['gene_id', 'gene_name']]
selected_sample = st.sidebar.selectbox("Select Sample Column", sample_cols)

# 3. Data Presentation
st.write(f"### Displaying counts for: {selected_sample}")
st.dataframe(df[['gene_id', 'gene_name', selected_sample]].head(20))

# 4. Visualization
st.write("### Gene Expression Distribution")
fig, ax = plt.subplots(figsize=(10, 5))
# Using your actual data columns to plot
sns.barplot(data=df.head(20), x='gene_name', y=selected_sample)
plt.xticks(rotation=90)
st.pyplot(fig)
