import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px

st.set_page_config(layout="wide")
st.title("🧬 DFU Biofilm: Multidimensional Interaction Assessment")

# 1. LOAD DATA
@st.cache_data
def load_data():
    return pd.read_csv("all_gene_counts_2.tsv", sep='\t')

df = load_data()

# 2. RESEARCHER INTERFACE
st.sidebar.header("Multidimensional Assessment")
gene_choice = st.sidebar.selectbox("Select Target Gene", df['gene_name'].unique())
selected_df = df[df['gene_name'] == gene_choice]

# 3. INTERACTION VISUALIZATION
col1, col2 = st.columns([1, 2])

with col1:
    st.subheader("Data Summary")
    st.write(f"Targeting: {gene_choice}")
    st.dataframe(selected_df.T)

with col2:
    st.subheader("Bacterial-Polymer Interaction Landscape")
    # Using Plotly for a professional, interactive 3D/Multivariate view
    fig = px.scatter_3d(selected_df.melt(id_vars=['gene_name']), 
                        x='variable', y='value', z='gene_name', 
                        color='value', title="Biofilm Response Dynamics")
    st.plotly_chart(fig, use_container_width=True)

# 4. STATISTICAL INTERACTION
st.subheader("Correlative Assessment")
fig2, ax2 = plt.subplots(figsize=(10, 4))
sns.heatmap(selected_df.select_dtypes(include=['float', 'int']).corr(), annot=True, cmap="viridis")
st.pyplot(fig2)
