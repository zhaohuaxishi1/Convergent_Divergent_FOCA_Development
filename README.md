# Convergent and divergent functional topographies of individualized human brain network and their developmental origins

E-mail: [zhaojianlong@mail.bnu.edu.cn](mailto:zhaojianlong@mail.bnu.edu.cn)

This repository provides code and source data that support the findings of the article entitled "*Convergent and divergent functional topographies of individualized human brain network and their developmental origins*“ 

## Oveview

This repository includes standalone software, source code, and demonstration data used in the analyses.
 All data required to reproduce the findings are publicly available, including:

- Individualized functional networks
- Group- and individual-level functional topographies
- Functional Topography Covariance (FOCA) matrices
- Source data for figure visualization

All files are hosted in a publicly accessible cloud repository (https://github.com/zhaohuaxishi1/Convergent_Divergent_FOCA_Development/tree/main).

### **Public datasets used**

- **MSC dataset:** https://openneuro.org/datasets/ds000224
- **dHCP dataset:** https://nda.nih.gov/
- **HCP dataset:** https://www.humanconnectome.org/

Source data supporting the findings are provided with this paper.

------

## Data

- **MSC_Group_Functional_Topography.dtseries.nii**
   Group-level functional topography from the MSC dataset, representing the canonical spatial organization of individualized functional networks.

  **HCP_Group_Functional_Topography.dtseries.nii**
   Group-level functional topography derived from the HCP adult dataset, used as a reference for adult network alignment and comparison.

  **dHCP_Group_Functional_Topography.dtseries.nii**
   Group-level functional topography derived from the dHCP neonatal dataset, representing the reference network architecture at birth.

  **MSC_FOCA_Hierarchical_Structure.mat**
   (Associated visualization: `./results/FOCA/MSC_FOCA_ReorderedMatrix.png`)
   Hierarchical organization of FOCA (Functional Topography Covariance) derived from the MSC dataset, reflecting nested antagonistic structures across cortical systems.

  **corrMatrix_dHCP_reordered.mat**
   (Associated visualization: `./results/FOCA/FOCA_ReorderedMatrix_dHCP.png`)
   FOCA matrix from the dHCP dataset, reordered according to hierarchical clustering for visualization and cross-age comparison.

  **37_44_week_dHCP_group_topography/**
   Group-level functional topographies across 37–44 weeks PMA, illustrating the progressive strengthening of negative connectivity within auditory and action-mode (AMN) networks, marking the emergence of antagonistic network organization during late gestation.

## **Code and Pipelines**

All analyses were conducted using open-source software and publicly available toolboxes.

### **Preprocessing**

- HCP minimal preprocessing pipeline: https://github.com/Washington-University/HCPpipelines/releases/
- dHCP structural pipeline: https://github.com/BioMedIA/dhcp-structural-pipeline
- dHCP functional pipeline: https://git.fmrib.ox.ac.uk/seanf/dhcp-neonatal-fmri-pipeline/-/tree/master
- dHCP fMRI Surface: https://git.fmrib.ox.ac.uk/seanf/dhcp-neonatal-fmri-pipeline/-/blob/master/doc/surface.md 

#### Postprocessing

- SPM12 toolbox: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
- GRETNA toolbox v2.0.0: https://www.nitrc.org/projects/gretna/
- CIFTI-Matlab toolbox v2: https://github.com/Washington-University/cifti-matlab/
- MATLAB R2020b: https://www.mathworks.com/products/matlab.html

#### Analysis

- Template Matching v1.0: https://github.com/DCAN-Labs/compare_matrices_to_assign_networks
- Support Vector Regression (SVR) model: https://github.com/ZaixuCui/Pattern_Regression_Clean and LIBSVM (v3.25): https://www.csie.ntu.edu.tw/~cjlin/libsvm/
- Generalized Additive Models (GAM): implemented in R using the mgcv package (https://cran.r-project.org/web/packages/mgcv/index.html)
- R 4.0.3: https://www.r-project.org
- MSM (Multimodal Surface Matching): https://github.com/ecr05/MSM_HOCR/releases

#### **Visualization**

- Hierarchical FOCA structure (MSC): `results/FOCA/MSC_FOCA_ReorderedMatrix.png`
- Hierarchical FOCA structure (dHCP): `results/FOCA/FOCA_ReorderedMatrix_dHCP.png`

- Neonatal 37–44 weeks group-level topographies: `/37_44_week_dHCP_group_topography/`
- Connectome Workbench: https://www.humanconnectome.org/software/connectome-workbench-R 
-  ggplot2: https://ggplot2.tidyverse.org/
- Surface Visualization:
  - HCP Adults:
     `S1200.L.very_inflated_MSMAll.32k_fs_LR.surf.gii`,
     `S1200.R.very_inflated_MSMAll.32k_fs_LR.surf.gii`
     Very-inflated cortical surfaces in the MSMAll (multimodal surface matching) 32k fs_LR space, used for adult topography visualization.
     [BALSA reference scene](https://balsa.wustl.edu/sceneFile/7qP5m)
  - dHCP Neonates:
     `week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii`,
     `week-40_hemi-right_space-dhcpSym_dens-32k_vinflated.surf.gii`
     Symmetric 40-week dHCP surfaces (32k space) for neonatal topography visualization.
     [Download template](https://biomedic.doc.ic.ac.uk/brain-development/downloads/dhcpSym_template.zip)

## **Custom Scripts and Analysis Steps**

Customized MATLAB and R scripts are provided for individualized network generation, FOCA computation, and visualization.
 Main analyses were performed following these steps:

1. **Generate group and individual functional networks:**
    `s1_1_generate_group_templates.m`, `s1_2_generate_individual_networks.m`
    Individualized networks generated using *Template Matching*.
2. **Compute individualized functional topography:**
    `s2_generate_individual_topography.m`
    Each network’s mean time series correlated with all cortical vertices.
3. **Construct FOCA matrix:**
    `s3_compute_foca_matrix.m`
    Quantifies convergent (positive) and divergent (negative) inter-network covariance.
4. **Assess reproducibility and variability:**
    `s4_consistency_variability_analysis.m`
    Quantifies cross-subject and within-subject FOCA stability.
5. **Hierarchical clustering of network organization:**
    `s5_hierarchical_clustering_and_modular_summary.m`
6. **Relate FOCA to cortical organizational features:**
    `s6_brain_axes_prediction_FOCA.m`
7. **Neonatal FOCA construction and comparison:**
    `s7_1_hierarchical_clustering_dHCP_comparison.m`,
    `s7_2_visualize_developmental_difference_FOCA_CohenD.m`
8. **Developmental modeling (GAM):**
    `s8_1_neonatal_FOCA_GAM_analysis.R`, `s8_2_age_effects_neonatal_FOCA_GAM.R`
9. **Predicting neurodevelopmental outcomes:**
    `s8_3_behavior_prediction_neonatal_FOCA_GLM.m`
    SVR models predicted 18-month cognitive and language outcomes.
10. **Comparing FOCA and conventional FC:**
     `s9_1_FOCA_vs_FC_stability_comparison.m`,
     `s9_2_FOCA_vs_FC_individual_dissimilarity.m`,
     `s9_3_NeonatalAgePrediction_SVR_FOCA_dHCP.m`,
     `s9_4_NeonatalAgePrediction_SVR_FC_dHCP.m`
     Compared stability, individual variability, and developmental sensitivity between FOCA and FC measures (with/without GSR).

------

✅ **Final Note:**
 This repository aims to facilitate reproducibility and future methodological extensions for functional topography and developmental connectomics research.
