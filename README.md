# Spatial Transcriptomics Analysis Platform

A web-based platform for analyzing and visualizing spatial transcriptomics data, specifically designed for 10x Visium data.

## Features

- **Data Upload & Processing**
  - Support for 10x Visium Spaceranger output
  - Direct .h5ad file upload
  - Pre-loaded example datasets (breast cancer, colon cancer)

- **Pre-processing**
  - Quality control metrics
  - Data normalization
  - Feature selection

- **Spatial Visualization**
  - Interactive spatial plots
  - Gene expression visualization
  - Customizable plotting parameters

- **Spatial Variable Features**
  - Moran's I statistics calculation
  - Highly variable gene identification
  - Spatial correlation analysis

## Data Preparation

### Processing Spaceranger Output
1. After running Spaceranger, locate your output directory
2. Run the preprocessing script:
   ```bash
   python process_SpaceRanger_out.py /path/to/spaceranger/output --output processed_data.h5ad
   ```
3. The generated `processed_data.h5ad` file can be uploaded directly through the web interface

## Installation

### Using Docker (Recommended)

1. Clone the repository: 