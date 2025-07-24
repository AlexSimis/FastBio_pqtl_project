# This script is used to jointly analyze the steady-state cis-pQTLs (TensorQTL) from all contexts.
# The main goal here is to (a) define the level of sharing between contexts and (b) identify cis-pQTLs with differential effects between contexs.

# Mashr (multivariate adaptive shrinkage)
library(mashr)
library(tidyverse)
library(data.table)
library(flashier)
path <- ""

# Full summary statistics.
pv1.stats <- fread(paste0(path, "PA1/proteins.1209/pa1.cis.sum.stats.txt"), header = T)
pv2.stats <- fread(paste0(path, "PA2/proteins.1209/pa2.cis.sum.stats.txt"), header = T)
nv1.stats <- fread(paste0(path, "NA1/proteins.1209/na1.cis.sum.stats.txt"), header = T)
nv2.stats <- fread(paste0(path, "NA2/proteins.1209/na2.cis.sum.stats.txt"), header = T)

# Top signals from the purmuted output of TensorQTL.
pv1.signals <- fread(paste0(path, "PA1/proteins.1209/pa1.cis_qtl.txt.gz"), header = T) 
pv2.signals <- fread(paste0(path, "PA2/proteins.1209/pa2.cis_qtl.txt.gz"), header = T)
nv1.signals <- fread(paste0(path, "NA1/proteins.1209/na1.cis_qtl.txt.gz"), header = T)
nv2.signals <- fread(paste0(path, "NA2/proteins.1209/na2.cis_qtl.txt.gz"), header = T)

###################################################
# Mash fit with strong (signals) and random samples.
# This approach is based on Matthew Stephens's eQTL 
# outline (https://stephenslab.github.io/mashr/articles/eQTL_outline.html).
###################################################

# IMPORTANT. For the downstream analysis, I have to (a)
# decide how to select the strong signals and (b) choose the 
# random subset of tests that will be used to estimate the correlation structure.

###################################################
# Preparing the input data.
###################################################

# Defining the matrix of Bhat, Shat per context.
pv1.stats |> select(phenotype_id, variant_id, slope, slope_se) |>
  mutate(gene_snp_pair = paste0(phenotype_id, "_", variant_id)) |>
  select(gene_snp_pair, slope, slope_se) -> pv1

pv2.stats |> select(phenotype_id, variant_id, slope, slope_se) |>
  mutate(gene_snp_pair = paste0(phenotype_id, "_", variant_id)) |>
  select(gene_snp_pair, slope, slope_se) -> pv2

nv1.stats |> select(phenotype_id, variant_id, slope, slope_se) |>
  mutate(gene_snp_pair = paste0(phenotype_id, "_", variant_id)) |>
  select(gene_snp_pair, slope, slope_se) -> nv1

nv2.stats |> select(phenotype_id, variant_id, slope, slope_se) |>
  mutate(gene_snp_pair = paste0(phenotype_id, "_", variant_id)) |>
  select(gene_snp_pair, slope, slope_se) -> nv2 

# Betas.
betas <- merge(select(pv1, gene_snp_pair, slope), select(pv2, gene_snp_pair, slope),
               by = "gene_snp_pair")
betas <- merge(betas, select(nv1, gene_snp_pair, slope),
               by = "gene_snp_pair")
betas <- merge(betas, select(nv2, gene_snp_pair, slope),
               by = "gene_snp_pair")

colnames(betas) <- c("gene_snp_pair", "b_pv1", "b_pv2", "b_nv1", "b_nv2")

# SEs.
se <- merge(select(pv1, gene_snp_pair, slope_se), select(pv2, gene_snp_pair, slope_se),
            by = "gene_snp_pair")
se <- merge(se, select(nv1, gene_snp_pair, slope_se),
            by = "gene_snp_pair")
se <- merge(se, select(nv2, gene_snp_pair, slope_se),
            by = "gene_snp_pair")

colnames(se) <- c("gene_snp_pair", "se_pv1", "se_pv2", "se_nv1", "se_nv2")

# Form the bhat matrix.
x.bhat <- betas[,-1]
x.bhat <- data.matrix(x.bhat)
rownames(x.bhat) <- betas$gene_snp_pair

# Form the shat matrix.
x.shat <- se[,-1]
x.shat <- data.matrix(x.shat)
rownames(x.shat) <- se$gene_snp_pair

# Save workspace.
# save.image("/bhat.shat.mash.RData")

# Setup the data for mash.
mash.data <- mash_set_data(x.bhat, x.shat)

###################################################
# (a) Picking strong signals.
###################################################

# Derive strong signals, across all contexts.
# Select the top SNP per Gene across all contexts.
bind_rows(pv1.stats, pv2.stats, nv1.stats, nv2.stats) |>
  mutate(group = rep(c("PV1", "PV2", "NV1", "NV2"),
                     times = c(nrow(pv1.stats), nrow(pv2.stats), nrow(nv1.stats), nrow(nv2.stats)))) |>
  group_by(phenotype_id) |>
  slice(which.min(pval_nominal)) -> pval.derived

top_pqtls <- paste(pval.derived$phenotype_id, pval.derived$variant_id, sep = "_")
rm(pval.derived); gc()

save.image("mash.model.Rdata")

###################################################
# (b) Selecting a random subset of tests.
###################################################
# Test the stability of Vhat between different random sumsamples.
# Random selection of tests.
Vhat.random.tests.list = data.frame()
Vhat = c()

vec = c(20000, 40000, 100000, 200000, 500000, 1304812, 2174686, nrow(mash.data$Bhat))
for (i in vec) {
  set.seed(12345)
  
  random.subset = sample(1:nrow(mash.data$Bhat), i)
  data.temp = mash_set_data(mash.data$Bhat[random.subset,], mash.data$Shat[random.subset,])
  Vhat = estimate_null_correlation_simple(data.temp)
  coord = data.frame(cor.PR = Vhat[2], cor.NR = Vhat[12])
  Vhat.random.tests.list = rbind(Vhat.random.tests.list, coord)
}

vec

Vhat.random.tests.list$subset = c("20k", "40k", "100k", "200k", "500k", "30%of.all.tests", "50%of.all.tests", "all.tests", "GTEx")
png("/Vhat.stability.png", width = 720, height = 720, res = 100);
plot(Vhat.random.tests.list$cor.PR, Vhat.random.tests.list$cor.NR,
     lines(Vhat.random.tests.list$cor.PR, Vhat.random.tests.list$cor.NR, col = "gray"),
     main = "Vhat stability of random subsamples", xlab = "Correlation PR1, PR2", ylab = "Correlation NR1, NR2");
text(Vhat.random.tests.list$cor.PR,Vhat.random.tests.list$cor.NR, Vhat.random.tests.list$subset, pos=1)
dev.off()

# From the Vhat stability plot I can see that the correlations cluster together after 500k randomly selected tests.

# Estimate the correlation structure from the random subset - null tests.
# Define the random subset of tests.
random.subset = 1:nrow(mash.data$Bhat) # Using all available tests.
data.temp = mash_set_data(mash.data$Bhat[random.subset,], mash.data$Shat[random.subset,])
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

Vhat
rm(pv1.stats, pv2.stats, nv1.stats, nv2.stats, x.bhat, x.shat, pv1, pv2, nv1, nv2)

# Setting up the data for the mash model.
data.random = mash_set_data(mash.data$Bhat[random.subset,], mash.data$Shat[random.subset,], V=Vhat)
data.strong = mash_set_data(mash.data$Bhat[top_pqtls,], mash.data$Shat[top_pqtls,], V=Vhat)

# Data-driven covariances from the strong signals.
U.pca = cov_pca(data.strong, 4) # PCA to estimate canditate covariance matrices.
U.ed = cov_ed(data.strong, U.pca) # Extreme deconvolution to refine the above estimates.

U.f = cov_flash(data.strong)

U.ed = cov_ed(data.strong, c(U.f, U.pca))

# Fit mash model on the random subset.
U.c = cov_canonical(data.random)
V.em = mash_estimate_corr_em(data.random, U.c, details = TRUE)

m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1) # Model fitment!

# Get the posterior estimates from the strong subset.
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

# Retrieve response pQTLs per group, i.e. pQTLs that are sig.
# in at least one timepoint and have a magnnitude of change greater than a threshold.
# My approach:
call.response.pqtls <- function(cols, pm, lfsr, group, thresh, mag) {
  
  lfsr.temp <- lfsr[, cols]
  pm.temp <- pm[, cols]
  colnames(pm.temp) <- c("pm_T1", "pm_T2")
  colnames(lfsr.temp) <- c("lfsr_T1", "lfsr_T2")
  
  print(identical(rownames(lfsr.temp), rownames(pm.temp)))
  
  cbind(lfsr.temp, pm.temp) |> data.frame() |> filter(lfsr_T1 < thresh | lfsr_T2 < thresh) |> mutate(sharing = ifelse((pm_T2/pm_T1) > (1/mag) & (pm_T2/pm_T1) < mag, "yes", "no"), fold.change = pm_T2/pm_T1) |> filter(sharing == "no") -> resp.pqtls
  
  assign(paste0(group, ".response.qtls"), resp.pqtls, envir = .GlobalEnv)
  
}

# Run full mashr model - compute posterior summaries - call response QTLs.
mash.model <- function(bhat, shat, cols, contexts, seed.number, strong.method, out.prefix){
  mash.data <- mash_set_data(bhat[,contexts], shat[,contexts])
  
  # Select 500k random tests - selected from the Vhat stability plot.
  #set.seed(seed.number)
  #random.subset = sample(1:nrow(mash.data$Bhat), 500000)
  random.subset = 1:nrow(mash.data$Bhat) # Using all the information of tests.
  # Derive strong signals, across all contexts.
  if (is.null(strong.method)){
    cat("Deriving strong signals based on nominal p-values ...")
    bind_rows(pv1.stats, pv2.stats, nv1.stats, nv2.stats) |>
      mutate(group = rep(c("PV1", "PV2", "NV1", "NV2"),
                         times = c(nrow(pv1.stats), nrow(pv2.stats), nrow(nv1.stats), nrow(nv2.stats)))) |>
      filter(group %in% cols) |>
      group_by(phenotype_id) |>
      slice(which.min(pval_nominal))-> pval.derived
    # Top signals.
    temp = pval.derived
  } else if (strong.method == "permuted.pval"){
    cat("Deriving strong signals based on the permuted p-values ...")
    
    # With this approach I'm using the information of permutations - hence I select signals based on the smallest permuted pvalue.
    bind_rows(pv1.signals, pv2.signals, nv1.signals, nv2.signals) |>
      mutate(group = rep(c("PV1", "PV2", "NV1", "NV2"),
                         times = c(nrow(pv1.signals), nrow(pv2.signals), nrow(nv1.signals), nrow(nv2.signals)))) |>
      group_by(phenotype_id) |>
      slice(which.min(pval_perm)) -> pval.perm.derived 
    # Top signals.
    temp = pval.perm.derived
  }
  
  # Check how the top pQTLs are distributed across contexts.
  print(table(temp$group))
  
  top_pqtls <- paste(temp$phenotype_id, temp$variant_id, sep = "_")
  rm(temp)
  
  # Estimate the correlation structure from the random subset.
  Vhat = estimate_null_correlation_simple( mash_set_data(mash.data$Bhat[random.subset,], mash.data$Shat[random.subset,]))
  cat("Vhat estimation based on the random subset", "\n"); 
  print(Vhat)
  write.table(Vhat, paste("vhat", out.prefix, "txt", sep = "."))
  
  # Setting up the data for the mash model.
  data.random = mash_set_data(mash.data$Bhat[random.subset,], mash.data$Shat[random.subset,], V=Vhat)
  data.strong = mash_set_data(mash.data$Bhat[top_pqtls,], mash.data$Shat[top_pqtls,], V=Vhat)
  
  # Data-driven covariances from the strong signals.
  U.pca = cov_pca(data.strong, 4) # PCA to estimate canditate covariance matrices.
  U.f = cov_flash(data.strong) # Use flash algorithm.
  U.f.ed = cov_ed(data.strong, c(U.f, U.pca))
  
  U.c = cov_canonical(data.random)
  
  # Fit mash model on the random subset.
  m = mash(data.random, Ulist = c(U.f.ed,U.c), outputlevel = 1) # Model fitment with flash!
  
  # Get the posterior estimates for the strong subset.
  m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
  
  pm <- get_pm(m2)
  lfsr <- get_lfsr(m2)
  
  # For PR individuals.
  call.response.pqtls(cols = c("b_pv1", "b_pv2"), pm = pm, lfsr = lfsr, group = "PR", thresh = 0.05, mag = 3)
  
  # For NR individuals.
  call.response.pqtls(cols = c("b_nv1", "b_nv2"), pm = pm, lfsr = lfsr, group = "NR", thresh = 0.05, mag = 1.5)
  
  # Store all response pqlts.
  resp.qtls <- c(rownames(PR.response.qtls), rownames(NR.response.qtls))
  resp.qtls <- unique(resp.qtls)
  
  # Check the lfsr of dr.qtls across groups.
  df <- data.frame(cbind(lfsr, pm))
  df <- df[resp.qtls, ]
  df$gene.snp <- rownames(df)
  
  colnames(df) <-  c("lfsrPA1", "lfsrPA2", "lfsrNA1", "lfsrNA2", 
                     "betaPA1", "betaPA2", "betaNA1", "betaNA2", "gene.snp.pair")
  
  df$discovered.in <- ifelse(df$gene.snp.pair%in%rownames(PR.response.qtls) & !df$gene.snp.pair%in%row.names(NR.response.qtls), "PR", 
                             ifelse(df$gene.snp.pair%in%rownames(NR.response.qtls) & !df$gene.snp.pair%in%rownames(PR.response.qtls), "NR", "shared"))
  
  
  #assign(paste("random.subset", out.prefix, sep = "."), random.subset, , envir = .GlobalEnv)
  #assign(paste("mash.data", out.prefix, sep = "."), mash.data, , envir = .GlobalEnv)
  #assign(paste("top.pqtls", out.prefix, sep = "."), top_pqtls, , envir = .GlobalEnv)
  assign(paste("mash.model", out.prefix, sep = "."), m2, , envir = .GlobalEnv)
  assign(paste("resp.pqtls", out.prefix, sep = "."), df, , envir = .GlobalEnv)
}

# Ground truth fitting mash with all tests!
mash.model(bhat = x.bhat, shat = x.shat, cols = c("PV1", "PV2", "NV1", "NV2"), strong.method = NULL, contexts = c(1:4),
           seed.number = 12345, out.prefix = "all.tests")  
# Save the workspace.
saveRDS(mash.model.all.tests, file = "/mnt/DataArray1/users/simistiras/pQTL_analysis/TensorQTL_pipeline/CxC_matched_samples/Results/Post_hoc/proteins.1209/mash.model.all.tests.rds")

# Load fitted model.
mash.model.all.tests <- readRDS("/mnt/DataArray1/users/simistiras/pQTL_analysis/TensorQTL_pipeline/CxC_matched_samples/Results/Post_hoc/proteins.1209/mash.model.all.tests.rds")

# Get posterior summaries.
pm <- get_pm(mash.model.all.tests)
lfsr <- get_lfsr(mash.model.all.tests)

# Get the number of significant conditions.
number.of.sig.pqtls  <- get_n_significant_conditions(mash.model.all.tests, thresh = 0.05)

number.of.sig.pqtls <- data.frame(gene.snp.pair = names(number.of.sig.pqtls), contexts = data.table(number.of.sig.pqtls)$number.of.sig.pqtls)

number.of.sig.pqtls |> filter(contexts > 0) -> number.of.sig.pqtls

colnames(pm) <- gsub("b_", "pm_", colnames(pm))

colnames(lfsr) <- gsub("b_", "lfsr_", colnames(lfsr))

bind_cols(pm, lfsr) |> mutate(gene.snp.pair = rownames(pm)) -> post.stats

# Attach the posterior stats.
mashr.cis.pqtls <- merge(post.stats, number.of.sig.pqtls, by = "gene.snp.pair")

write.table(mashr.cis.pqtls, "/mashr.cis.pqtls.txt", row.names = F, sep = "\t")

# Call response (or differential) cis-pQTLs between time points for each dietary groups.
# For PR individuals.
call.response.pqtls(cols = c("b_pv1", "b_pv2"), pm = pm, lfsr = lfsr, group = "PR", thresh = 0.05, mag = 1.5)

# For NR individuals.
call.response.pqtls(cols = c("b_nv1", "b_nv2"), pm = pm, lfsr = lfsr, group = "NR", thresh = 0.05, mag = 1.5)

# Store all response pqlts.
resp.qtls <- c(rownames(PR.response.qtls), rownames(NR.response.qtls))
resp.qtls <- unique(resp.qtls)

# Check the lfsr of dr.qtls across groups.
df <- data.frame(cbind(lfsr, pm))
#df <- data.frame(pm = m2$result$PosteriorMean, lfsr = m2$result$lfsr)
df <- df[resp.qtls, ]
df$gene.snp <- rownames(df)

colnames(df) <-  c("lfsrPA1", "lfsrPA2", "lfsrNA1", "lfsrNA2", 
                   "betaPA1", "betaPA2", "betaNA1", "betaNA2", "gene.snp.pair")

# Define in which context the cis-pQTLs were discovered.
df$discovered.in <- ifelse(df$gene.snp.pair%in%rownames(PR.response.qtls) & !df$gene.snp.pair%in%row.names(NR.response.qtls), "PR", 
                           ifelse(df$gene.snp.pair%in%rownames(NR.response.qtls) & !df$gene.snp.pair%in%rownames(PR.response.qtls), "NR", "shared"))

