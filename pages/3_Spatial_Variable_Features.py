import streamlit as st
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import squidpy as sq
import libpysal
from esda.moran import Moran

st.set_page_config(
    page_title="Spatial Variable Features",
)

st.title("Spatial Variable Features")

# Get adata from session state
if 'adata' not in st.session_state:
    st.error('Please upload data in the Home page!')
    st.stop()

adata = st.session_state.adata

# Check if Moran's I has been calculated
has_moran = 'moran_i' in adata.var.columns

if not has_moran:
    # Ask whether to calculate Moran's I
    if st.button("Calculate Moran's I for spatial variable features"):
        with st.spinner("Calculating Moran's I statistics..."):
            try:
                # Create spatial weights matrix
                coords = adata.obsm['spatial']
                kdt = libpysal.weights.KNN.from_array(coords, k=8)
                w = kdt.weights
                
                # Calculate Moran's I for each gene
                moran_scores = []
                moran_pvals = []
                
                for gene in adata.var_names:
                    gene_exp = adata[:, gene].X.flatten()
                    mi = Moran(gene_exp, kdt)
                    moran_scores.append(mi.I)
                    moran_pvals.append(mi.p_sim)
                
                # Add results to adata.var
                adata.var['moran_i'] = moran_scores
                adata.var['moran_i_pval'] = moran_pvals
                
                has_moran = True
                st.success("Moran's I calculation completed!")
                
                # Show some statistics
                st.write("Top 10 spatially variable genes by Moran's I:")
                top_spatial_genes = adata.var.sort_values('moran_i', ascending=False).head(10)
                st.write(top_spatial_genes[['moran_i', 'moran_i_pval']])
                
            except Exception as e:
                st.error(f"An error occurred: {str(e)}")

# If there is Moran's I information, show download button
if has_moran:
    if st.download_button(
        label="Download processed AnnData object",
        data=adata.write_h5ad(),
        file_name="spatial_variable_features.h5ad",
        mime="application/x-hdf5"
    ):
        st.success("Download started!")

# Add progress bar
with st.spinner('Calculating spatial variable features...'):
    # Use 'seurat' method instead of 'seurat_v3'
    sc.pp.highly_variable_genes(adata, n_top_genes=1000, flavor='seurat')
    
    # Get variable gene results
    var_genes = pd.DataFrame({
        'Gene': adata.var_names,
        'Variability': adata.var['dispersions_norm'],
        'highly_variable': adata.var['highly_variable']
    })
    
    # Sort by variability
    var_genes = var_genes.sort_values('Variability', ascending=False)

# Show top 10 results
st.subheader('Top 10 Variable Genes')
st.write(
    var_genes[['Gene', 'Variability']]
    .head(10)
)

# Visualize variability of top 10 genes
fig, ax = plt.subplots(figsize=(10, 6))
sns.barplot(data=var_genes.head(10), x='Gene', y='Variability')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
st.pyplot(fig)


    