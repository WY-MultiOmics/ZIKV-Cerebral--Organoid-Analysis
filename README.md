# Zika Virus-Induced HNRNP Lactylation Impairs Neural Development in Human Brain Organoids

This repository contains the full multi-omics data analysis pipeline, source code, and visualization scripts for the manuscript:

> **Yu Wang#, Tao Chen#, Wen Zhao, Zheyu Xu, Yuanyuan Zhang, Yue Xiong, Zhujiang, Yunhuang Yang, Maili Liu, Pu Chen*, Rui Hu*. "Zika Virus-Induced HNRNP Lactylation Impairs Neural Development in Human Brain Organoids". (2026)**
> 
> *# Co-first authors; \* Co-corresponding authors (hurui@apm.ac.cn, puchen@whu.edu.cn)*

---
## Project Overview
This repository contains the analytical pipelines used to investigate the macromolecular underpinnings of Zika virus (ZIKV) infection in human brain organoids (hBOs). 

Our research systematically maps the molecular kinetic cascades triggered by ZIKV, revealing that **RNA splicing** is the most sensitive vulnerability during neurogenesis. We identify the **heterogeneous nuclear ribonucleoprotein (hnRNP) family** as the central regulatory hub, and demonstrate that **site-specific lysine lactylation** of HNRNPD and HNRNPU is essential for maintaining their RNA-chaperoning functions and structural integrity during human cortical organogenesis.

## Key Research Methodologies
- **Multi-Omics Integration**: Time-resolved proteomics, Lactylome profiling, and single-cell RNA sequencing (scRNA-seq).
- **Structural & Biochemical Validation**: Hydrogen-Deuterium Exchange Mass Spectrometry (HDX-MS), NMR titration, and CRISPR/site-directed mutagenesis.
- **Model Systems**: Human brain organoids (hBOs).

## Pipeline Organization
- `01_Time_resolved_Proteomics/`: Computational framework for temporal clustering and functional enrichment of ZIKV-induced proteomic changes.
- `02_lactylation/`: Analysis of the lactylome landscape and site-specific modification dynamics of the hnRNP family.
- `03_single_cell/`: End-to-end scRNA-seq workflow including integration, cell-type annotation, and trajectory inference (RNA velocity).

## 📁 Repository Structure

The analytical workflow is structured chronologically according to the multi-omics layers described in the paper:

```text
my-paper-code/
├── README.md
├── environment.yml
└── scripts/
    ├── 01_proteomics/
    │   ├── 01_data_preprocessing.R
    │   ├── 02_temporal_clustering.R
    │   ├── 03_time_functional_go.R
    │   ├── 04_time_functional_kegg.R
    │   ├── 05_time_point_dep_volcano_curve.R
    │   ├── 06_dep_kegg_down_plot.R
    │   ├── 07_dep_go_up_plots.R
    │   ├── 08_ferroptosis_zscore_heatmaps.R
    │   ├── 09_ferroptosis_mechanism_profiles.R
    │   ├── 10_mock_timepoint_mfuzz.R
    │   ├── 11_mock_cluster_go_enrichment.R
    │   ├── 12_mock_cluster_time_heatmaps.R
    │   └── 13_zikv_cluster_time_heatmaps.R
    ├── 02_lactylation/
    │   ├── 01_lactylome_dep_volcano.R
    │   ├── 02_global_dep_volcano.R
    │   ├── 03_kla_site_dep_heatmap.R
    │   ├── 04_kegg_up_enrichment.R
    │   ├── 05_hnrnp_kla_sites_heatmap.R
    │   ├── 06_hnrnp_family_kla_intensify.R
    │   ├── 07_go_up_enrichment.R
    │   └── 08_kla_extract_21mer_motif.R
    └── 03_single_cell/
        ├── 01_mock1_qc.ipynb
        ├── 02_per_sample_qc.ipynb
        ├── 03_scvi_integration.ipynb
        ├── 04_clustering_and_markers.ipynb
        ├── 05_cell_type_analysis.ipynb
        ├── 06_celltypist_refinement.ipynb
        ├── 07_zikv_infected_mapping.ipynb
        ├── 08_infected_cell_annotations.ipynb
        ├── 09_infected_cell_analysis.ipynb
        ├── 11_zikv_infected_deg_analysis.ipynb
        ├── 12_deg_go_analysis.ipynb
        ├── 13_zikv_mock_rna_velocity_combined.ipynb
        ├── 14_loom_filtering_alignment.ipynb
        ├── 14_scvelo_analysis.ipynb
        ├── 15_rna_velocity_visualization.ipynb
        ├── 16_zikv_mock_trajectory_pipeline.ipynb
        ├── 17_zikv_mock_group_trajectory.ipynb
        ├── 18_velocity_magnitude_statistics.ipynb
        └── 19_multi_gene_velocity_dynamics.ipynb
        
##  Requirements & Environment
- **Proteomics & Lactylation:** R (v4.3.0+)
The Single cell analyses were performed in a Linux environment. The software dependencies are managed using Conda.
- **Single-cell Analysis:** Python (v3.10+) with `scanpy`, `scvelo`, `scvi-tools`, and `pyscenic`.
To replicate the analysis environment, run:
```bash
conda env create -f environment_r.yml

##  Data Availability
scRNA-seq: Raw sequencing data is deposited at GEO (Accession: GSE338563).
Mass Spectrometry: The time-resolved and lactylation proteomics datasets have been deposited in the ProteomeXchange Consortium via the PRIDE partner repository under accession number (PXD081059, PXD080910)

## Contact
For questions regarding the code or pipeline, please contact:Rui Hu (hurui@apm.ac.cn)