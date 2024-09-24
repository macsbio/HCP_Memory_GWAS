## Beyond single variants: A pathway-based approach to explore the genetic basis of memory

This repository contains all the code used to analyze the GWAS data of the human connectome project to study genetic variants related to memory performance. 
The publication is currently under review.

Content:
* ```1_remove_duplicates.R```: Contains the code needed to remove duplicate loci that are present on both chips from the Human Connectome Project while keeping the duplicated locus with the highest call rate
* ```2_quality_control_preprocess.sh```: Standard quality control of the genotype data; calculate genomic PCA and genetic relationship matrix
* ```3_create_phenotype_files.R```: Takes the phenotype data from the Human Connectome Project and filters it; Also calculates the test correlations and the residuals
* ```4_gwas.sh```: The code to run the GWAS using GCTA
* ```5_gwas2FUMA.R```: Getting the output from GCTA into a format readable by FUMA
* ```5a_qq_manhattan.R```: For plotting Manhattan and QQ plots.
* ```6_FUMA2Networks.R```: Taking the FUMA output, filtering it, and visualizing it in Cytoscape networks
