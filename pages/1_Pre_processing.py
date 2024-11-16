import streamlit as st
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pickle
import io
import os
import tempfile
import squidpy as sq
import pandas as pd
from io import BytesIO

# Check if 'leiden' exists in adata.obs columns
def check_preprocessing(adata):
    if 'leiden' not in adata.obs.columns:
        return False
    return True

# Preprocessing function
def run_preprocessing(adata):
    # Make gene names unique before saving raw data
    adata.var_names_make_unique()
    
    # Save raw data
    adata.raw = adata
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)
    
    return adata

# Main preprocessing interface
if 'adata' in st.session_state:
    adata = st.session_state.adata
    
    if not check_preprocessing(adata):
        st.warning("Data needs preprocessing")
        if st.button("Start Preprocessing"):
            with st.spinner("Preprocessing in progress..."):
                adata = run_preprocessing(adata)
                st.session_state.adata = adata
            st.success("Preprocessing completed!")
            
            # Create temporary file and save data
            with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp:
                adata.write_h5ad(tmp.name)
                with open(tmp.name, 'rb') as f:
                    bytes_data = f.read()
                os.unlink(tmp.name)  # Delete the temporary file
                
            if st.download_button(
                label="Save Preprocessed Data",
                data=bytes_data,
                file_name="preprocessed_data.h5ad",
                mime="application/x-h5"
            ):
                st.success("Data saved successfully!")
            
            # Create tabs for different visualizations
            tab1, tab2, tab3 = st.tabs(["UMAP Clustering Results", "Spatial Transcriptomics Clustering", "Feature Gene Analysis"])
            
            with tab1:
                fig1, ax1 = plt.subplots(figsize=(10, 8))
                sc.pl.umap(adata, color=['leiden'], show=False, ax=ax1)
                st.pyplot(fig1)
                
            with tab2:
                # Create two columns for side-by-side plots
                col1, col2 = st.columns(2)
                
                with col1:
                    # Plot H&E image only
                    fig_he, ax_he = plt.subplots(figsize=(8, 8))
                    sq.pl.spatial_scatter(adata, color="leiden", img=True, ax=ax_he, alpha=0)
                    st.pyplot(fig_he)
                    
                with col2:
                    # Plot spatial transcriptomics with leiden clusters
                    fig_spatial, ax_spatial = plt.subplots(figsize=(8, 8))
                    sq.pl.spatial_scatter(adata, color="leiden", img=False, ax=ax_spatial)
                    st.pyplot(fig_spatial)
                
            with tab3:
                # Check if rank_genes_groups has already been calculated
                if 'rank_genes_groups' not in adata.uns:
                    with st.spinner("Calculating feature genes..."):
                        sc.tl.rank_genes_groups(
                            adata,
                            groupby='leiden',
                            method='wilcoxon',
                            n_genes=10,
                            use_raw=True,
                            tie_correct=True
                        )
                        # Update session state after calculation
                        st.session_state.adata = adata
                
                # 将结果转换为DataFrame格式
                result = adata.uns['rank_genes_groups']
                groups = result['names'].dtype.names
                
                # 创建一个字典来存储每个cluster的基因
                cluster_genes = {}
                for group in groups:
                    genes = result['names'][group][:10]  # 获取前10个基因
                    scores = result['scores'][group][:10]  # 获取对应的分数
                    pvals = result['pvals'][group][:10]   # 获取对应的p值
                    
                    # 创建该cluster的DataFrame
                    cluster_df = pd.DataFrame({
                        'Gene': genes,
                        'Score': scores,
                        'P-value': pvals
                    })
                    
                    # 在Streamlit中显示每个cluster的结果
                    st.subheader(f'Cluster {group} Top 10 Marker Genes')
                    st.dataframe(cluster_df.set_index('Gene').reset_index())
    else:
        st.success("Data has been preprocessed")
        
        # Create temporary file and save data
        with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp:
            adata.write_h5ad(tmp.name)
            with open(tmp.name, 'rb') as f:
                bytes_data = f.read()
            os.unlink(tmp.name)  # Delete the temporary file
            
        if st.download_button(
            label="Save Preprocessed Data",
            data=bytes_data,
            file_name="preprocessed_data.h5ad",
            mime="application/x-h5"
        ):
            st.success("Data saved successfully!")
        
        # Create tabs for different visualizations
        tab1, tab2, tab3 = st.tabs(["UMAP Clustering Results", "Spatial Transcriptomics Clustering", "Feature Gene Analysis"])
        
        with tab1:
            fig1, ax1 = plt.subplots(figsize=(10, 8))
            sc.pl.umap(adata, color=['leiden'], show=False, ax=ax1)
            st.pyplot(fig1)
            
        with tab2:
            # Create two columns for side-by-side plots
            col1, col2 = st.columns(2)
            
            with col1:
                # Plot H&E image only
                fig_he, ax_he = plt.subplots(figsize=(8, 8))
                sq.pl.spatial_scatter(adata, color="leiden", img=True, ax=ax_he, alpha=0)
                st.pyplot(fig_he)
                
            with col2:
                # Plot spatial transcriptomics with leiden clusters
                fig_spatial, ax_spatial = plt.subplots(figsize=(8, 8))
                sq.pl.spatial_scatter(adata, color="leiden", img=False, ax=ax_spatial)
                st.pyplot(fig_spatial)
            
            with tab3:
                # Check if rank_genes_groups has already been calculated
                if 'rank_genes_groups' not in adata.uns:
                    with st.spinner("Calculating feature genes..."):
                        sc.tl.rank_genes_groups(
                            adata,
                            groupby='leiden',
                            method='wilcoxon',
                            n_genes=10,
                            use_raw=True,
                            tie_correct=True
                        )
                        # Update session state after calculation
                        st.session_state.adata = adata
                
                # 将结果转换为DataFrame格式
                result = adata.uns['rank_genes_groups']
                groups = result['names'].dtype.names
                
                # 创建一个字典来存储每个cluster的基因
                cluster_genes = {}
                for group in groups:
                    genes = result['names'][group][:10]  # 获取前10个基因
                    scores = result['scores'][group][:10]  # 获取对应的分数
                    pvals = result['pvals'][group][:10]   # 获取对应的p值
                    
                    # 创建该cluster的DataFrame
                    cluster_df = pd.DataFrame({
                        'Gene': genes,
                        'Score': scores,
                        'P-value': pvals
                    })
                    
                    # 在Streamlit中显示每个cluster的结果
                    st.subheader(f'Cluster {group} Top 10 Marker Genes')
                    st.dataframe(cluster_df.set_index('Gene').reset_index())
else:
    st.error("Please upload a data file first")


