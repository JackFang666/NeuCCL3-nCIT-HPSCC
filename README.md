# Analysis Code for: CCL3+ Neutrophil Signature Predicts Response to nCIT in HPSCC

This repository contains the R code and analysis pipeline used in the manuscript:
**"CCL3+ neutrophil signature predicts response to neoadjuvant toripalimab plus chemotherapy in hypopharyngeal squamous cell carcinoma: A phase II trial"**

## 📂 Repository Structure
* `Code/`: Contains the R scripts for single-cell and bulk RNA-seq analysis.
* `Data/`: Contains processed metadata and small input files (Raw data is hosted externally).

## 🧬 Data Availability
The raw single-cell RNA-seq and bulk RNA-seq data generated in this study have been deposited in the **GSA-Human** database.
* Note: Please download the raw matrices from GSA/GEO and place them in the `./data/raw` folder to run the scripts.

## 💻 System Requirements
* **OS**: Windows 10/11 or Linux / macOS
* **Software**: R version >= 4.1.0
* **Key Dependencies**: 
  * Seurat (v4.0 or v5.0)
  * Harmony
  * CellChat
  * clusterProfiler

## 🚀 How to Run
1. Clone this repository.
2. Download raw data from [Link to Data].
3. Run `HPSCC_Analysis_Pipeline.R` to reproduce the clustering and signature analysis.

