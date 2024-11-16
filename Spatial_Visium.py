import streamlit as st
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import squidpy as sq
import tempfile

# Add dropdown menu for example datasets
example_datasets = {
    "None": None,
    "Breast Cancer": "Examples/breastcancer.h5ad",
    "Colon Cancer": "Examples/coloncancer.h5ad"
}

col1, col2 = st.columns(2)

with col1:
    # Add file uploader for h5ad files
    uploaded_file = st.file_uploader("Upload a .h5ad file:", type=['h5ad'])

with col2:
    selected_example = st.selectbox(
        "Or select an example dataset:",
        options=list(example_datasets.keys())
    )
    load_example = st.button("Load Example", disabled=(selected_example == "None"))

# Handle data loading
if uploaded_file is not None:
    # Read uploaded h5ad file
    with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp:
        tmp.write(uploaded_file.getvalue())
        adata = sc.read_h5ad(tmp.name)
        st.session_state['adata'] = adata
    os.unlink(tmp.name)  # Clean up temp file
elif load_example and selected_example != "None":
    # Load example dataset
    example_path = example_datasets[selected_example]
    if os.path.exists(example_path):
        adata = sc.read_h5ad(example_path)
        st.session_state['adata'] = adata
        st.success(f"Successfully loaded {selected_example} dataset!")
    else:
        st.error(f"Example file not found: {example_path}")



