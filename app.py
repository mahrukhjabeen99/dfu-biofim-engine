import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

st.title("🧬 DFU Biofilm Research Engine")

try:
    # Use 'sep=None' and 'engine=python' to automatically detect if it's comma or tab separated
    df = pd.read_csv("authentic_results.csv", sep=None, engine='python')
    
    st.write("### Quantitative Gene Expression Results")
    st.dataframe(df.head(20)) # Show first 20 rows

    # Check if 'Log2FC' exists
    if 'Log2FC' in df.columns:
        st.write("### Gene Expression Fold Change Distribution")
        fig, ax = plt.subplots(figsize=(10, 5))
        sns.histplot(df['Log2FC'].dropna(), bins=30, kde=True, color='purple')
        st.pyplot(fig)
    else:
        st.error(f"The column 'Log2FC' was not found. Available columns: {list(df.columns)}")
        
except Exception as e:
    st.error(f"Error loading CSV: {e}")
