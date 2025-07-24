# This is a shell script to run cis-qtl mapping using tensorQTL.
# The script produces (1) gene-level pQTLs, (2) full summary statistics, and (3) independent cis-pQTLs per gene.

# Import parameters:
path=$1
bfile=$2
phenotype=$3
outprefix=$4
covs=$5

# With gene-level permutations.
python3 -m tensorqtl $path/$bfile \
    $path/$phenotype \
    $path/$outprefix \
    --covariates $path/$covs --seed 123456 --mode cis

# Full summary statistics.
python3 -m tensorqtl $path/$bfile \
    $path/$phenotype \
    $path/$outprefix \
    --covariates $path/$covs --seed 123456 --mode cis_nominal

# Independent cis-pQTLs per gene.
cis_output=".cis_qtl.txt.gz"
python3 -m tensorqtl $path/$bfile \
    $path/$phenotype \
    $path/$outprefix \
    --covariates $path/$covs \
    --cis_output $path/$outprefix$cis_output --seed 123456 --mode cis_independent
