import streamlit as st
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import squidpy as sq
import requests
import json
import tempfile
import os

st.set_page_config(
    page_title="Spatial Transcriptomics Visualization",
    layout="wide"
)

st.title("Spatial Transcriptomics Data Visualization")

def query_gene_function(gene_name):
    """Query gene function using Ollama API"""
    try:
        response = requests.post(
            'http://100.108.16.21:11434/api/generate',
            json={
                'model': 'medllama2',
                'prompt': f"What is the function of the {gene_name} gene? Please provide a brief description."
            },
            timeout=10,
            stream=False
        )
        if response.status_code == 200:
            # Split response by lines and process each line as separate JSON
            response_text = ''
            for line in response.text.strip().split('\n'):
                try:
                    if line.strip():  # Skip empty lines
                        json_response = json.loads(line)
                        if 'response' in json_response:
                            response_text += json_response['response']
                except json.JSONDecodeError:
                    continue
            return response_text if response_text else 'No information available'
        return f"Error: Unable to fetch gene information (Status code: {response.status_code})"
    except Exception as e:
        return f"Error connecting to Ollama service: {str(e)}"

if 'adata' not in st.session_state:
    st.warning('Please upload data file on the home page first!')
else:
    adata = st.session_state.adata
    
    # Create five-column layout instead of four
    col1, col2, col3, col4, col5 = st.columns([0.5, 1, 1, 1, 1])
    
    with col1:
        st.subheader("Settings")
        
        # Gene input
        gene_name = st.text_input("Input Gene of Interest")
        
        # Basic visualization settings
        color_map = 'RdYlBu_r'
        spot_size = st.slider(
            "Set Spot Size",
            min_value=150,
            max_value=300,
            value=200
        )
        
        if gene_name:
            if adata.raw is not None:
                available_genes = adata.raw.var_names
                adata_use = adata.raw
            else:
                available_genes = adata.var_names
                adata_use = adata
            
            if gene_name not in available_genes:
                st.error(f"Gene '{gene_name}' not found in dataset, please check the gene name.")
            else:
                # Get gene expression max value
                gene_max = float(adata_use[:, gene_name].X.max())
                
                # Color range slider
                color_range = st.slider(
                    "Set Color Range",
                    min_value=0.0,
                    max_value=gene_max,
                    value=(0.0, gene_max),
                    step=float(gene_max/100),
                    key=f"color_range_{gene_name}"
                )

    # Move all visualization code outside of col1
    if gene_name and gene_name in available_genes:
        # Visualization columns (col2 to col5)
        with col2:
            st.subheader("Gene Expression Violin Plot")
            fig1, ax1 = plt.subplots(figsize=(8, 6))
            sc.pl.violin(
                adata,
                gene_name,
                groupby='leiden',
                ax=ax1,
                use_raw=True
            )
            st.pyplot(fig1)
            plt.close(fig1)

        with col3:
            st.subheader("Spatial Gene Expression")
            fig2, ax2 = plt.subplots(figsize=(6, 6))
            sc.pl.spatial(
                adata,
                color=gene_name,
                spot_size=spot_size,
                color_map=color_map,
                ax=ax2,
                use_raw=True,
                img_key=None,
                vmin=color_range[0],
                vmax=color_range[1]
            )
            st.pyplot(fig2)
            plt.close(fig2)

        with col4:
            st.subheader("Spatial Distribution by Leiden")
            fig3, ax3 = plt.subplots(figsize=(6, 6))
            sc.pl.spatial(
                adata,
                color='leiden',
                spot_size=spot_size,
                ax=ax3,
                img_key=None
            )
            st.pyplot(fig3)
            plt.close(fig3)

        with col5:
            st.subheader("H&E Background Image")
            fig4, ax4 = plt.subplots(figsize=(6, 6))
            sq.pl.spatial_scatter(adata, color="leiden", img=True, ax=ax4, alpha=0)
            st.pyplot(fig4)
            plt.close(fig4)

        # Add gene function query after all columns
        st.write("### Gene Function (MedLLaMA2 output):")
        with st.spinner('Querying gene function...'):
            gene_function = query_gene_function(gene_name)
            st.write(gene_function)

        # Update session state
        st.session_state.prev_gene_name = gene_name
        st.session_state.prev_spot_size = spot_size