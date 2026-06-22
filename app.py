import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

st.title("🧬 DFU Biofilm Research Engine")

# 1. LOAD YOUR REAL DATA
@st.cache_data
def load_data():
    return pd.read_csv("all_gene_counts_2.tsv", sep='\t')

df = load_data()

# 2. RESEARCHER INPUTS
st.sidebar.header("Filter Results")
# Assuming your file has columns for 'Bacteria' and 'Polymer'
# If these columns don't exist, we must use your Metadata file to map them
selected_bacteria = st.sidebar.selectbox("Bacterial Species", df['Bacteria'].unique())
selected_polymer = st.sidebar.selectbox("Polymer Treatment", df['Polymer'].unique())

# 3. AUTHENTIC DATA RETRIEVAL
# We filter the database for the exact experimental condition
results = df[(df['Bacteria'] == selected_bacteria) & (df['Polymer'] == selected_polymer)]

st.write(f"### Observed Transcriptomic Counts")
st.dataframe(results)

# 4. SCIENTIFIC VISUALIZATION
st.write("### Gene Expression Response")
fig, ax = plt.subplots()
# We plot the raw count values as measured in your lab
sns.barplot(data=results.head(15), x='gene_id', y='counts')
plt.xticks(rotation=90)
st.pyplot(fig)
