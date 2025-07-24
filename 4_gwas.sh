#!bin/bash
./gcta/gcta64 --mlma-loco \
--bfile merged_plink_qc \
--grm merged_plink_qc_grm \
--pheno pheno_iwrd_residuals.txt \
--out gwas_output/iwrd \
--thread-num 10

./gcta/gcta64 --mlma-loco \
--bfile merged_plink_qc \
--grm merged_plink_qc_grm \
--pheno pheno_pics_residuals.txt \
--out gwas_output/pics \
--thread-num 10

./gcta/gcta64 --mlma-loco \
--bfile merged_plink_qc \
--grm merged_plink_qc_grm \
--pheno pheno_list_residuals.txt \
--out gwas_output/list \
--thread-num 10
