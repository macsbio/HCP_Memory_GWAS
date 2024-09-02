#!bin/bash
plink1.9 --bfile merged_plink --geno 0.05 --mind 0.05 --maf 0.05 --hwe 1e-6 --make-bed --out merged_plink_qc
plink1.9 --bfile merged_plink_qc --pca --out merged_plink_qc_pca
./gcta/gcta64 --bfile merged_plink_qc --make-grm --thread-num 10 --out merged_plink_qc_grm
./gcta/gcta64 --grm merged_plink_qc_grm --make-bK-sparse 0.05 --out merged_plink_qc_spgrm
