import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import h5py
import json
import os
import argparse
from datetime import datetime

def process_visium_data(input_dir, output_file=None):
    """
    Process 10x Visium data from spaceranger output
    
    Parameters:
    -----------
    input_dir : str
        Path to spaceranger output directory
    output_file : str, optional
        Path to save the processed data file. If None, will save in the input directory
    
    Returns:
    --------
    adata : AnnData object
        Processed data object
    """
    
    # Check if input directory exists
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"Input directory {input_dir} does not exist!")
    
    # Set output file
    if output_file is None:
        sample_name = os.path.basename(input_dir)
        output_file = os.path.join(input_dir, f"{sample_name}_processed.h5ad")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir:  # 只在输出路径包含目录时创建目录
        os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Read the Visium data
        print(f"Reading data from {input_dir}...")
        adata = sc.read_visium(input_dir)
        
        # Add some basic information
        adata.uns['processed_date'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Save the processed data
        print(f"Saving processed data to {output_file}...")
        adata.write_h5ad(output_file)
        
        # Print basic information about the dataset
        print("\nData processing completed successfully!")
        print(f"Number of spots: {adata.n_obs}")
        print(f"Number of genes: {adata.n_vars}")
        print(f"Output file: {output_file}")
        
        return adata
        
    except Exception as e:
        print(f"Error occurred during processing: {str(e)}")
        raise

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process 10x Visium data from spaceranger output')
    parser.add_argument('input_dir', type=str, help='Path to spaceranger output directory')
    parser.add_argument('--output', type=str, help='Path to save the processed data file', default=None)
    
    # Parse arguments
    args = parser.parse_args()
    
    # Process the data
    adata = process_visium_data(args.input_dir, args.output) 