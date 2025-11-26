# Computational Environment for Figure Generation

All figures were generated using R version 4.3.1 on macOS (arm64 architecture). Figure generation relied on tidyverse packages (ggplot2, dplyr, tidyr, readr), Bioconductor genomic tools (GenomicRanges, BSgenome, IRanges), and specialized visualization libraries (ggrepel, ggbreak, pheatmap, viridis). The exact package versions and full R session information are provided below to ensure full computational reproducibility.

R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS 26.0

## Software Environment Used for Figure Generation

The figures in this repository were generated using the following R packages and versions.

### Core plotting and data handling
- ggplot2 3.5.2  
- ggrepel 0.9.6  
- ggnewscale 0.5.1  
- ggbreak 0.1.4  
- pheatmap 1.0.12  
- reshape2 1.4.4  
- readxl 1.4.5  
- dplyr 1.1.4  
- tidyr 1.3.1  
- stringr 1.5.1  
- data.table 1.17.2  
- viridis 0.6.5  
- viridisLite 0.4.2  
- openxlsx 4.2.8  
- readr 2.1.5  

### Genomic context (used for locus annotation and genomic visualizations)
- GenomicRanges 1.52.1  
- GenomeInfoDb 1.36.4  
- IRanges 2.34.1  
- S4Vectors 0.38.2  
- BiocGenerics 0.46.0  
- regioneR 1.32.0  
- Biostrings 2.68.1  
- GenomicAlignments 1.36.0  
- rtracklayer 1.60.1  

### Supporting utilities
- forcats 1.0.0  
- conflicted 1.2.0  
- yaml 2.3.10  
- scales 1.4.0  
- patchwork 1.3.0  
- Matrix 1.6-5  
- MatrixGenerics 1.12.3  
- SummarizedExperiment 1.30.2  
- XVector 0.40.0  
- Rsamtools 2.16.0  
- tibble 3.2.1  
- purrr 1.0.4  
- memoise 2.0.1  
- hms 1.1.3  

*Namespace-only dependencies such as `rlang`, `vctrs`, `pillar`, and others were loaded automatically and are part of standard tidyverse/Bioconductor infrastructure.*


