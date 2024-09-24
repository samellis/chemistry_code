"""
streamlit_ketcher_app
A skeleton app for drawing molecules with streamlit_ketcher.

To run:
streamlit run streamlit_ketcher_app.py
Sam Ellis
"""

import streamlit as st
from streamlit_ketcher import st_ketcher

# Page layout
st.set_page_config(layout='wide', page_title="Molecule Searcher with Streamlit")
st.title("Molecule Drawing with Streamlit")

st.text("Draw a molecule and press apply to display the SMILES.")
smiles = st_ketcher()

if smiles:
    # Whenever st_ketcher updates
    st.text(f"Molecule drawn is {smiles}.")