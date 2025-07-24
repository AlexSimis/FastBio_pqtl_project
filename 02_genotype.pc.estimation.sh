#!/bin/sh
# Genotype PCs estimation to control for population structure. The estimation of genotype PCs will be based on a subset of independent SNPs.

# PR individuals.
plink_v2.00a --bfile /mnt/DataArray1/users/simistiras/pQTL_analysis/TensorQTL_pipeline/CxC_matched_samples/Data/Genotypes/PA1/PA1.genotypes --keep pa1.matched.ids.txt --make-bed --out PA1.genotypes.matched

# Prune SNPs to estimate the genotype PCs.
plink_v2.00a --bfile /mnt/DataArray1/users/simistiras/pQTL_analysis/TensorQTL_pipeline/CxC_matched_samples/Data/Genotypes/PA1/PA1.genotypes.matched --indep-pairwise 250 50 0.2 --bad-ld --out indep.PA1.genotypes.matched

# Estimate 3 geno PCs.
plink_v2.00a --bfile /mnt/DataArray1/users/simistiras/pQTL_analysis/TensorQTL_pipeline/CxC_matched_samples/Data/Genotypes/PA1/PA1.genotypes.matched --extract indep.PA1.genotypes.matched.prune.in --pca 3 --out /mnt/DataArray1/users/simistiras/pQTL_analysis/TensorQTL_pipeline/CxC_matched_samples/Data/Genotypes/PA1/3geno.pcs.PA1

# NR individuals.
# Prune SNPs to estimate the genotype PCs.
plink_v2.00a --bfile /mnt/DataArray1/users/simistiras/pQTL_analysis/TensorQTL_pipeline/CxC_matched_samples/Data/Genotypes/NA1/NA1.genotypes --keep na1.matched.ids.txt --make-bed --out NA1.genotypes.matched

plink_v2.00a --bfile /mnt/DataArray1/users/simistiras/pQTL_analysis/TensorQTL_pipeline/CxC_matched_samples/Data/Genotypes/NA1/NA1.genotypes.matched --indep-pairwise 250 50 0.2 --bad-ld --out indep.NA1.genotypes.matched

# Estimate 3 geno PCs.
plink_v2.00a --bfile /mnt/DataArray1/users/simistiras/pQTL_analysis/TensorQTL_pipeline/CxC_matched_samples/Data/Genotypes/NA1/NA1.genotypes.matched --extract indep.NA1.genotypes.matched.prune.in --pca 3 --out /mnt/DataArray1/users/simistiras/pQTL_analysis/TensorQTL_pipeline/CxC_matched_samples/Data/Genotypes/NA1/indep.3geno.pcs.NA1
