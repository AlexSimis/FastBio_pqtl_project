# This script is used to assess the replication of the FastBio top cis-pQTLs against the UKB-PPP.
# Cis-summary statistics of UKB-PPP have downloaded through "https://www.synapse.org/Synapse:syn51365303" EUR discovery.

library(tidyverse)
library(data.table)

path <- "../../CxC_mapping/Results/"
# Load the FastBio top cis-pQTLs.
pa1.stats <- fread(paste0(path, "pa1.cis_qtl.txt.gz"))
pa2.stats <- fread(paste0(path, "pa2.cis_qtl.txt.gz"))
na1.stats <- fread(paste0(path, "na1.cis_qtl.txt.gz"))
na2.stats <- fread(paste0(path, "na2.cis_qtl.txt.gz"))

# Store them in 1 df.
bind_rows(pa1.stats, pa2.stats, na1.stats, na2.stats) %>%
  mutate(context = rep(c("PR1", "PR2", "NR1", "NR2"), times = c(nrow(pa1.stats), nrow(pa2.stats), nrow(na1.stats), nrow(na2.stats))), gene.snp.pair = paste(phenotype_id, variant_id, sep = "_")) %>% 
  filter(qval < 0.05) -> cis.pqtls

# Remove obsolete dfs.
rm(pa1.stats, pa2.stats, na1.stats, na2.stats)
gc()

# Load cis-pQTLs from UKB-PPP pre-filtered at 1e-5.
pqtls.ukb <- fread("../Post_hoc_analysis/replication/data/pqtls.ukbb.1e5.filtered.txt.gz")

# Sanity checks.
glimpse(pqtls.ukb)

table(pqtls.ukb$CHROM)

pqtls.ukb <- pqtls.ukb[!which(pqtls.ukb$CHROM == "CHROM"),]

colnames(pqtls.ukb)[10] <- "gene"

# Set proper data types.
temp <- apply(pqtls.ukb[, -c(3,4,10)], 2,
              function(x) as.numeric(x))

pqtls.ukb <- cbind(pqtls.ukb[,c(3,4,10)], temp)

# Setting variables.
pqtls.ukb %>% 
  mutate(pval = 10^(-LOG10P), variant_id = paste(CHROM, GENPOS, ALLELE0,
                                                 ALLELE1, sep = ":")) %>%
  filter(pval < 0.05) %>% 
  mutate(gene.snp.pair = paste(gene, variant_id, sep = "_")) -> pqtls.ukb

nrow(distinct(pqtls.ukb, gene)) # Sanity check.

rm(temp)

# Due to an error while retrieving the ukb summary stats the second part of the double id genes is missing.
# Track the proper names of missing genes.
to.replace <- setdiff(indep.pqtls$phenotype_id, pqtls.ukb$gene)

pqtls.ukb %>% 
  filter(pqtls.ukb$gene %in% gsub("_.*", "", to.replace)) %>% 
  group_by(gene) %>% 
  summarise(n())

rep(to.replace[c(4,3,5,2,1)], each = c(60846,2009,1596,1978,696	))

pqtls.ukb[which(pqtls.ukb$gene %in% gsub("_.*", "", to.replace)), "gene"] <- rep(to.replace[c(4,3,5,2,1)], times = c(60846,2009,1596,1978,696	))

pqtls.ukb %>% 
  filter(pqtls.ukb$gene %in% to.replace) %>% 
  group_by(gene) %>% 
  summarise(n()) # Sanity check.

pqtls.ukb %>% 
  mutate(temp = variant_id) %>% 
  separate(temp, into = c("chr", "pos")) %>% 
  mutate(chr.pos = paste(chr, pos, sep = ":")) -> pqtls.ukb

pqtls.ukb %>% 
  distinct(gene)

# Merge cis-pqtls against ukb-pqtls at 1e-5.
cis.pqtls %>% 
  group_by(context) %>% 
  summarise(n())

cis.pqtls %>% 
  group_by(phenotype_id) %>% 
  dplyr::slice(which.min(pval_perm))

# Merge FastBio pQTLs with UKB-sum stats based on the same gene-snp pair.
merged.cis.pqtls <- merge(cis.pqtls, pqtls.ukb, by = "gene.snp.pair")

merged.cis.pqtls %>% 
  group_by(context) %>% 
  filter(BETA*slope < 0)

plot(merged.cis.pqtls$BETA, merged.cis.pqtls$slope)
plot(merged.cis.pqtls$A1FREQ, merged.cis.pqtls$af)

# Define a Bonferroni corrected replication threshold.

# Keep the ones with concordant direction.
merged.cis.pqtls %>% 
  filter(BETA*slope > 0) -> merged.cis.pqtls

# Replication rate of all top cis-pQTLs.
merged.cis.pqtls %>%
  filter(pval < (0.05/630)) %>% 
  distinct(phenotype_id) %>% 
  summarise(n()/630)

merged.cis.pqtls %>% 
  group_by(phenotype_id) %>% 
  dplyr::slice(which.min(pval_perm)) -> plot.data

# Check the pearson's cor. coef. of betas between studies.
cor(plot.data$BETA, plot.data$slope)

