#!bin/bash
./gcta/gcta64 --bfile merged_plink_qc --grm-sparse merged_plink_qc_spgrm --fastGWA-mlm --pheno pheno_iwrd_residuals.txt --thread-num 10 --out gene_assoc_iwrd_residuals

./gcta/gcta64 --bfile merged_plink_qc --grm-sparse merged_plink_qc_spgrm --fastGWA-mlm --pheno pheno_pics_residuals.txt --thread-num 10 --out gene_assoc_pics_residuals

./gcta/gcta64 --bfile merged_plink_qc --grm-sparse merged_plink_qc_spgrm --fastGWA-mlm --pheno pheno_list_residuals.txt --thread-num 10 --out gene_assoc_list_residuals
