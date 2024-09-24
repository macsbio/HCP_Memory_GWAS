This repository contains all the code used in the study "Beyond single variants: A pathway-based approach to explore the genetic basis of memory" 2024.

Content:
    File 1      Contains the code needed to remove duplicate loci that are present on both chips from the Human
                Connectome Project while keeping the duplicated locus with the highest call rate

    File 2      Standard quality control of the genotype data
                Calculate genomic PCA and genetic relationship matrix
    
    File 3      Takes the phenotype data from the Human Connectome Project and filters it
                Also calculates the test correlations and the residuals
    
    File 4      The code to run the GWAS using GCTA
    
    File 5      Getting the output from GCTA into a format readable by FUMA
                File 5a for plotting Manhattan and QQ plots.
    
    File 6      Taking the FUMA output, filtering it, and visualizing it in Cytoscape networks
