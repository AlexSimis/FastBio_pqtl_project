# This script is used to generate and select covariates to adjust in cis-pQTL mapping.
library(tidyverse)
library(caroline)
library(data.table)

# Loading data.
args <- commandArgs(trailingOnly = TRUE)

qced.phenotype.data <- args[1] # Import QCed phenotype data.
meta.data <- args[2]

qced.data <- fread(qced.phenotype.data)

# Load the metadata.
meta <- readr::read_delim(meta.data,
                          col_names = TRUE,
                          show_col_types = FALSE)
colnames(meta)[1] = 'SampleID' # Set proper column name.

# Optional additional metadata.
# meta <- fread("../../../Biomarkers_GWAS/data/FastBio_biomarkers_full_data.csv")
# colnames(meta)[1] = 'SampleID' # Set proper column name.

# Define covariates.
meta$age2 <- meta$age^2
meta$sex <- as.factor(meta$sex)
meta$fasting <- as.factor(meta$fasting)
meta$timePoint <- as.factor(meta$timePoint)
meta$medication_use <- as.factor(meta$medication_use)
meta$Smoking <- as.factor(meta$Smoking)
meta$Alcohol <- as.factor(meta$Alcohol)
glimpse(meta)

# Make a unique SampleID, remember I am aggregating the observations from timepoints.
paired.ids <- separate(select(qced.data,SampleID), SampleID, into = c('ids', 'timePoint'), sep = '_')
paired.ids$SampleID <- qced.data$SampleID

paired.ids %>% 
  filter(duplicated(ids)) -> paired.ids # Participants who attended both time points.

qced.data.paired <- qced.data[grepl(paste(paired.ids$ids, collapse = "|"), qced.data$SampleID),]

table(qced.data.paired$fasting, qced.data.paired$timePoint) # Sanity check.

# Isolate SampleIDs per Diet - Timepoint group.
pa1.ids <- qced.data.paired[which(qced.data.paired$fasting == "F" & qced.data.paired$timePoint == "T1"),1]

pa2.ids <- qced.data.paired[which(qced.data.paired$fasting == "F" & qced.data.paired$timePoint == "T2"),1]

na1.ids <- qced.data.paired[which(qced.data.paired$fasting == "NF" & qced.data.paired$timePoint == "T1"),1]

na2.ids <- qced.data.paired[which(qced.data.paired$fasting == "NF" & qced.data.paired$timePoint == "T2"),1]

# write.table(cbind(pa1.ids, pa1.ids), "Data/Matching.ids/pa1.matched.ids.txt", quote = F, row.names = F, col.names = F, sep = "\t")
# write.table(cbind(pa2.ids, pa2.ids), "Data/Matching.ids/pa2.matched.ids.txt", quote = F, row.names = F, col.names = F, sep = "\t")
# write.table(cbind(na1.ids, na1.ids), "Data/Matching.ids/na1.matched.ids.txt", quote = F, row.names = F, col.names = F, sep = "\t")
# write.table(cbind(na2.ids, na2.ids), "Data/Matching.ids/na2.matched.ids.txt", quote = F, row.names = F, col.names = F, sep = "\t")

rm(qced.data); gc()

# Make a covariate file containing both protein and genotype derived PCs.
# Load the genotype PCs (calculated via plink) to add them as covariate.
to.load <- list.files("Data/Covariates/",
                      pattern = "indep.3geno.pcs.*.eigenvec",
                      full.names = T)
geno.pcs <- lapply(to.load, function(x){
  fread(x) %>% 
    select(!c(`#FID`))
})

# List names.
names(geno.pcs)[1] <- "NR1"
names(geno.pcs)[2] <- "NR2"
names(geno.pcs)[3] <- "PR1"
names(geno.pcs)[4] <- "PR2"

# But since we have targeted QC on phenotypes and genotypes we want to keep only 
# the intersect of the SampleIDs - individuals.
# Check the missingness.
matching.pheno.geno <- function(pheno, geno.pcs, out.prefix){
  to.keep <- intersect(pheno$SampleID, geno.pcs$IID)
  
  # Drop the missing participant from the qced data.
  pheno <- filter(pheno, SampleID%in%to.keep)
  nrow(pheno) # sanity check.
  
  # Check the arrangement and the sampleID match
  geno.pcs <-  geno.pcs[match(pheno$SampleID, geno.pcs$IID),]
  colnames(geno.pcs)[1] = "SampleID"
  
  assign(paste0("pheno.", out.prefix), pheno, envir = .GlobalEnv)
  assign(paste0("geno.pcs.", out.prefix), geno.pcs, envir = .GlobalEnv)
  
}

matching.pheno.geno(pheno = filter(qced.data.paired, SampleID%in%pa1.ids$SampleID), 
                    geno.pcs = geno.pcs$PR1, 
                    out.prefix = "PA1")

matching.pheno.geno(pheno = filter(qced.data.paired, SampleID%in%pa2.ids$SampleID), 
                    geno.pcs = geno.pcs$PR2, 
                    out.prefix = "PA2")

matching.pheno.geno(pheno = filter(qced.data.paired, SampleID%in%na1.ids$SampleID), 
                    geno.pcs = geno.pcs$NR1, 
                    out.prefix = "NA1")

matching.pheno.geno(pheno = filter(qced.data.paired, SampleID%in%na2.ids$SampleID), 
                    geno.pcs = geno.pcs$NR2, 
                    out.prefix = "NA2")

# Data manipulation to meet the needs of PCAtools.
# Metadata : SampleIDs as rows.
meta.data <- select(meta, age, age2, GC, BMI, medianIsize, SBP, DBP, sex, Smoking,
                    Alcohol, medication_use)
meta.data %>% 
  mutate(across(everything(), as.numeric)) -> meta.data

rownames(meta.data) <- meta$SampleID

# Make a function.
pca.format <- function(data, meta, out.prefix){
  library(PCAtools)
  meta <- meta[which(rownames(meta)%in%data$SampleID),]
  rownames(meta) <- data$SampleID
  t.data <- t(as.matrix(dplyr::select(data, starts_with(c('Sa','OID')))))
  table(is.na(t.data))
  t.data <- data.frame(t.data)
  colnames(t.data ) <- t.data [1,]
  t.data  <- t.data [-1,]
  t.data[,] <- sapply(t.data, as.numeric)
  p <- pca(t.data, metadata = meta, scale = T)
  out.prefix <-  assign(out.prefix, p, envir = .GlobalEnv)
}

# Perform PCA.
pca.format(data = pheno.PA1, meta = meta.data, out.prefix='pca.pa1')
pca.format(data = pheno.PA2, meta = meta.data, out.prefix='pca.pa2')
pca.format(data = pheno.NA1, meta = meta.data, out.prefix='pca.na1')
pca.format(data = pheno.NA2, meta = meta.data, out.prefix='pca.na2')

# Adjusting for the number of protein PCs as derived from the elbow method & 3 genotype PCs plus sex, age, BMI & medication use.
# Isolate the phenotype eigenvectors.
pa1.pcs <- pca.pa1$rotated[,1:findElbowPoint(pca.pa1$variance)]
pa2.pcs <- pca.pa2$rotated[,1:findElbowPoint(pca.pa2$variance)]
na1.pcs <- pca.na1$rotated[,1:findElbowPoint(pca.na1$variance)]
na2.pcs <- pca.na2$rotated[,1:findElbowPoint(pca.na2$variance)]

identical(rownames(pa1.pcs), geno.pcs.PA1$SampleID) # Sanity check.
identical(rownames(pa1.pcs), pheno.PA1$SampleID)

meta$sex <- as.numeric(meta$sex)
meta$medication_use <- as.numeric(meta$medication_use)

# Transpose & merge all covariates.
transpose.df <- function(data1, data2, data3, named.df) {
  data1$SampleID <- rownames(data1)
  covs <- merge(data1, data2, by = "SampleID")
  covs <- merge(covs, data3, by = "SampleID")
  
  temp <- t(as.matrix(covs))
  temp[,1:5]
  temp <- data.frame(temp)
  colnames(temp) <- temp[1,]
  temp <- temp[-1,]
  temp[,] <- sapply(temp, as.numeric)
  assign(named.df, temp, envir=.GlobalEnv)
}

# PR1
transpose.df(data1 = pa1.pcs, 
             data2 = select(meta, SampleID, age, sex, BMI, 
                            medication_use),
             data3 = geno.pcs.PA1, named.df = 'covs.df.pa1')
glimpse(covs.df.pa1)

# PR2
transpose.df(data1 = pa2.pcs, 
             data2 = select(meta, SampleID, age, sex, BMI, 
                            medication_use),
             data3 = geno.pcs.PA2, named.df = 'covs.df.pa2')
# NR1
transpose.df(data1 = na1.pcs, 
             data2 = select(meta, SampleID, age, sex, BMI,
                            medication_use),
             data3 = geno.pcs.NA1, named.df = 'covs.df.na1')
# NR2
transpose.df(data1 = na2.pcs, 
             data2 = select(meta, SampleID, age, sex, BMI, 
                            medication_use),
             data3 = geno.pcs.NA2, named.df = 'covs.df.na2')

# Export covariates per group.
path = ""
prefix = ""

write.table(covs.df.pa1, paste(paste0(path,'covariates.PA1'), 
                               length(colnames(covs.df.pa1)), 
                               nrow(covs.df.pa1), prefix, 'txt', sep = '.'), row.names = T, col.names = T, quote = F, sep = '\t')
write.table(covs.df.pa2, paste(paste0(path,'covariates.PA2'), 
                               length(colnames(covs.df.pa2)), 
                               nrow(covs.df.pa2), prefix, 'txt', sep = '.'), row.names = T, col.names = T, quote = F, sep = '\t')
write.table(covs.df.na1, paste(paste0(path,'covariates.NA1'), 
                               length(colnames(covs.df.na1)), 
                               nrow(covs.df.na1), prefix, 'txt', sep = '.'), row.names = T, col.names = T, quote = F, sep = '\t')
write.table(covs.df.na2, paste(paste0(path,'covariates.NA2'), 
                               length(colnames(covs.df.na2)), 
                               nrow(covs.df.na2), prefix, 'txt', sep = '.'), row.names = T, col.names = T, quote = F, sep = '\t')
