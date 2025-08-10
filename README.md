# SPCA (Spatial Pathway and Class Analysis)
Software tools for pathway and chemical class analysis for Mass spectrometry imaging data. The data analysis tools are built to be compatible with MetaSpace Datasets. There are currently 2 pathway and chemical class analysis implementations for MSI data: Kernel PCA and Z-score.
    - Kernel PCA info: https://doi.org/10.1186/s12859-022-05005-1
    - Z-score info: https://doi.org/10.1371/journal.pcbi.1000217

There are 3 tutorials: 
  - SPA_tutorial.ipynb: Process metaspace dataset, segment dataset, load and filter reactome, run Kernel PCA pathway analysis, cluster pathway scores and analyze intercluster variance
  - SPA_spatclust_tutorial.ipynb: Process metaspace dataset, segment dataset, load and filter reactome, run Z-score pathway analysis, analyze pathway scores using variability and Moran's i
  - SCA_tutorial.ipynb: Process metaspace dataset, segment dataset, run Kernel PCA Chemical class analysis, analyze class scores using variability and Moran's i
