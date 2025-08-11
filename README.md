# SPCA (Spatial Pathway and Class Analysis)
Software tools for pathway and chemical class analysis for Mass spectrometry imaging data. The data analysis tools are built to be compatible with MetaSpace annotated datasets. There are currently 2 pathway and chemical class analysis implementations for MSI data: Kernel PCA and Z-score.
    - Kernel PCA info: https://doi.org/10.1186/s12859-022-05005-1
    - Z-score info: https://doi.org/10.1371/journal.pcbi.1000217

To install, I recommend creating a virtual environment in order to keep the packages straight and not mess with the python system.
"""python3.9 -m venv myenv"""

Activate the environment
"""source myenv/bin/activate   # macOS/Linux"""
"""myenv\Scripts\activate      # Windows"""

Then install the dependencies (ensure you are in the same directory as the requirements.txt file:
"""pip install -r requirements.txt"""


There are 3 tutorials: 
  - SPA_tutorial.ipynb: Process metaspace dataset, segment dataset, load and filter reactome, run Kernel PCA pathway analysis, cluster pathway scores and analyze intercluster variance
  - SPA_spatclust_tutorial.ipynb: Process metaspace dataset, segment dataset, load and filter reactome, run Z-score pathway analysis, analyze pathway scores using variability and Moran's i
  - SCA_tutorial.ipynb: Process metaspace dataset, segment dataset, run Kernel PCA Chemical class analysis, analyze class scores using variability and Moran's i
