import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

st.title("🧬 DFU Biofilm Research Engine")
st.markdown("### Real-Time Transcriptomic Analysis of DFU Biofilm")

# Load your actual research data
@st.cache_data
def load_data():
    return pd.read_csv("authentic_results.csv")

df = load_data()

# Display Data
st.write("### Quantitative Gene Expression Results")
st.dataframe(df)

# Visualize Data
st.write("### Gene Expression Fold Change Distribution")
fig, ax = plt.subplots(figsize=(10, 5))
sns.histplot(df['Log2FC'], bins=30, kde=True, color='purple')
st.pyplot(fig)
