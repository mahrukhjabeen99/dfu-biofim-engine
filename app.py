
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

st.title("🔬 DFU Biofilm Analytics Engine")

# Load data
if os.path.exists("counts.csv"):
    df = pd.read_csv("counts.csv")
    st.success("Data loaded successfully.")
else:
    st.warning("Please upload counts.csv to GitHub.")

# Inputs
staph = st.number_input("S. aureus Load", value=1000)
zeta = st.slider("Zeta Potential (mV)", -80.0, 80.0, 0.0)

# Physics calculation
dist = np.linspace(1, 50, 100)
v_total = -1.5e-20 / (12 * np.pi * (dist * 1e-9)**2) + (52.0 * zeta * 0.15 * np.exp(-1.2 * dist))

# Plot
fig, ax = plt.subplots()
ax.plot(dist, v_total)
st.pyplot(fig)
