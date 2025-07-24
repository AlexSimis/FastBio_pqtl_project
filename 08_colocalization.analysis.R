# This script is used to perform (a) phewas queries and identify associated GWAS traits with 
# differential cis-QTLs and (b) based on these results perform colocalization in each context.

library(data.table)
library(coloc)
library(tidyverse)
library(ieugwasr)
path = ""
# Load response pQTLs.
resp.pqtls <- fread(paste0(path, "/pQTL_analysis/TensorQTL_pipeline/CxC_mapping/Post_hoc_analysis/mashr/results/response.pqtls.txt"))

# Gene name.
resp.pqtls %>%
  separate(gene.snp.pair, into = c("phenotype_id", "variant_id"), sep = "_") %>%
  select(phenotype_id, variant_id) -> resp.pqtls

# Load context-specific cis-pQTL on FGF21.
fgf21.pqtl <- fread("../../pQTL_analysis/TensorQTL_pipeline/CxC_mapping/Results/pa2.cis_qtl.txt.gz") %>% 
  filter(phenotype_id == "FGF21") %>% 
  select(phenotype_id, variant_id)

# Add the FGF21 associated cis-pQTL with the list of differential cis-pQTLs,
resp.pqtls <- rbind(resp.pqtls, fgf21.pqtl)

# Fetch proxies - (r2 > 0.8 using plink)
indep.pqtl.proxies <- fread(paste0(path, "/pQTL_analysis/TensorQTL_pipeline/CxC_mapping/Results/indep.snps.proxies.txt.ld"))
resp.pqtls <- merge(resp.pqtls, select(indep.pqtl.proxies, SNP_A, SNP_B), by.x = "variant_id", by.y = "SNP_A")

# Track the lead = diff. qtl and the snp = proxy.
colnames(resp.pqtls)[match(c("variant_id", "SNP_B"), colnames(resp.pqtls))] <- c("lead.snp", "snp")
rm(indep.pqtl.proxies, fgf21.pqtl); gc()

# Load TensorQTL summary stats.
pa1.stats <- fread(paste0(path, "/pQTL_analysis/TensorQTL_pipeline/CxC_mapping/Results/pa1.cis.sum.stats.txt"))
pa2.stats <- fread(paste0(path, "/pQTL_analysis/TensorQTL_pipeline/CxC_mapping/Results/pa2.cis.sum.stats.txt"))
na1.stats <- fread(paste0(path, "/pQTL_analysis/TensorQTL_pipeline/CxC_mapping/Results/na1.cis.sum.stats.txt"))
na2.stats <- fread(paste0(path, "/pQTL_analysis/TensorQTL_pipeline/CxC_mapping/Results/na2.cis.sum.stats.txt"))

# Store them in 1 df.
bind_rows(pa1.stats, pa2.stats, na1.stats, na2.stats) %>%
  mutate(context = rep(c("PR1", "PR2", "NR1", "NR2"), times = c(nrow(pa1.stats), nrow(pa2.stats), nrow(na1.stats), nrow(na2.stats))),
         gene.snp.pair = paste(phenotype_id, variant_id, sep = "_")) -> pqtl.sum.stats

# Select the pGenes of interest.
pqtl.sum.stats %>%
  filter(phenotype_id %in% resp.pqtls$phenotype_id) -> pqtl.sum.stats

# Remove obsolete dfs.
rm(pa1.stats, pa2.stats, na1.stats, na2.stats)
gc()

# Get rsIDs.
rsids <- fread("../Data/Genotypes/rsIDs_annotated_snps.txt")

# Add the rsIDs.
pqtl.sum.stats <- merge(pqtl.sum.stats, select(rsids, hg38, hg19, RefSNP_id),
                        by.x = c("variant_id"),
                        by.y = c("hg38"), all.x = T
)

pqtl.sum.stats %>% 
  rename(rsID = RefSNP_id) -> pqtl.sum.stats

rm(rsids); gc()

# List of candidate resp.pQTLs for PheWAS.
merge(resp.pqtls, select(pqtl.sum.stats, variant_id, hg19, rsID), by.x = "snp", by.y = "variant_id") %>%
  distinct(snp, .keep_all = T) %>%
  mutate(temp = hg19) %>%
  separate(temp, into = c("chr", "position", "nea", "ea")) -> resp.pqtls

# Load the annotated genes. Note the TSS that I have used is on hg38 but I need 
# to retrieve the TSS on hg19 because of "ieugwasr" package.
genes.tss <- fread("../Data/Phenotypes/Annotated.OlinkIDs.1500.csv")

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# I have to lift over the TSS to hg37.
mart.hg37 <- useDataset("hsapiens_gene_ensembl", useMart("ensembl", host = "https://grch37.ensembl.org"))
# Hg38 build.
mart.hg38 <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

listAttributes(mart.hg37)
geneTSS_hg37 <- getBM(
  attributes = c(
    "ensembl_gene_id", "ensembl_transcript_id",
    "hgnc_symbol", "chromosome_name",
    "transcription_start_site",
    "gene_biotype", "strand", "start_position",
    "end_position"
  ),
  values = TRUE, mart = mart.hg37
) # hg37

# Note in hg37 has not information of conical transcript.
geneTSS_hg38 <- getBM(
  attributes = c(
    "ensembl_gene_id", "ensembl_transcript_id",
    "hgnc_symbol", "chromosome_name",
    "transcript_is_canonical", "transcription_start_site",
    "gene_biotype", "strand"
  ),
  values = TRUE, mart = mart.hg38
) # hg38 "transcript_is_canonical"

# Sanity check.
geneTSS_hg38 %>%
  filter(hgnc_symbol %in% resp.pqtls$phenotype_id & transcript_is_canonical == 1) -> test

# geneTSS_hg38 %>%
#   filter(hgnc_symbol %in% genes.tss$Assay & transcript_is_canonical == 1) -> test

temp <- merge(test, genes.tss, by.x = "hgnc_symbol", by.y = "Assay")

plot(temp$transcription_start_site.x, temp$transcription_start_site.y) # Sanity check.

# Take the transcript that I know is canonical from the hg38.
geneTSS_hg37 %>%
  filter(ensembl_transcript_id %in% temp$ensembl_transcript_id) -> geneTSS_hg37

# Remove obsolete dfs.
rm(geneTSS_hg38, mart.hg37, mart.hg38, temp, test, genes.tss, ensembl)

detach("package:biomaRt", unload = TRUE)

# Run phewas.
# Select the response pGene with the associated summary stats from the context
phewas.fun <- function(data, group, gene, thresh) {
  # Chr pos columns.
  data %>%
    filter(context == group & phenotype_id == gene) %>%
    mutate(tmp = hg19) %>%
    separate(tmp, into = c("chr", "position"), sep = ":") -> temp
  
  # Attach the TSS of the gene in hg37.
  resp.pqtls <- merge(resp.pqtls, geneTSS_hg37, by.x = "phenotype_id", by.y = "hgnc_symbol", all.x = T)
  to.select <- c("rsID", "phenotype_id", "lead.snp", "transcription_start_site", "ea", "nea", "chr", "position")
  
  temp.phewas <- merge(temp, select(resp.pqtls, to.select), by = c("rsID", "phenotype_id"))
  
  temp.phewas$lead.snp <- temp.phewas[which(temp.phewas$variant_id == temp.phewas$lead.snp), "rsID"]
  
  # Perform phewas with default parameters.
  phewas.res <- phewas(variants = temp.phewas$rsID)
  
  # Curate traits.
  phewas.res %>%
    filter(!is.na(n)) %>% # Drop missing sample size.
    group_by(trait) %>% # Keep the most sig.association per trait
    slice(which.max(abs(beta / se))) -> phewas.res
  
  # Keep track of the lead snp & the queried SNP.
  phewas.res %>%
    mutate(
      lead.snp = ifelse(rsid == unique(temp.phewas$lead.snp), unique(temp.phewas$lead.snp),
                        unique(temp.phewas$lead.snp)
      ),
      queried.snp = rsid
    ) %>%
    filter(p < thresh) %>%
    arrange(p) -> phewas.res
  
  # Attach the pQTL data.
  to.select <- c("rsID", "phenotype_id", "transcription_start_site")
  merge(dplyr::select(temp.phewas, to.select), phewas.res, by.x = "rsID", by.y = "rsid", suffix = c(".pQTL", ".trait")) %>%
    mutate(
      region = paste0(chr, ":", transcription_start_site - 1000000 , "-", transcription_start_site + 1000000 ),
      snp.window = paste0(chr, ":", position - 1000000, "-", position + 1000000)
    ) %>% 
    arrange(p) %>% 
    filter(!grepl("eqtl", id)) -> phewas.res
  
  assign("pqtl.stats", temp, envir = .GlobalEnv)
  assign("phewas.res", phewas.res, envir = .GlobalEnv)
}

# Run the function.
phewas.fun(data = pqtl.sum.stats, group = "PR2", gene = "LBR", thresh = 1e-5)
phewas.fun(data = pqtl.sum.stats, group = "PR2", gene = "FGF21", thresh = 1e-5)
# Retrieve associations.
coloc.res <- lapply(1:nrow(phewas.res), function(i) {
  # Fetch associations
  # assoc.temp <- associations(variants = phewas.res$region[i], id = phewas.res$id[i])
  assoc.temp <- tryCatch(
    {
      associations(variants = phewas.res$region[i], id = phewas.res$id[i])
      # associations(variants = phewas.res$snp.window[i], id = phewas.res$id[i])
    },
    error = function(e) {
      message(paste("Error fetching associations for", phewas.res$snp.window[i], "-", conditionMessage(e)))
      return(NULL)
    },
    warning = function(w) {
      message(paste("Warning for", phewas.res$snp.window[i], "-", conditionMessage(w)))
      invokeRestart("muffleWarning") # suppress if you want
    }
  )
  
  # Skip to next if associations failed
  if (!is.data.frame(assoc.temp) || nrow(assoc.temp) == 0) {
    return(list(pqtl.assoc = NULL, coloc.res = NULL))
  }
  
  assoc.temp %>%
    filter(!is.na(beta)) %>%
    group_by(id, trait) %>%
    distinct(rsid, .keep_all = T) -> assoc.temp
  
  assoc.temp$queried.snp <- phewas.res$queried.snp[i]
  assoc.temp$gene <- phewas.res$phenotype_id[i]
  assoc.temp$lead.snp <- phewas.res$lead.snp[i]
  
  # Retrieve trait/study info
  trait.info <- tryCatch(
    {
      gwasinfo(unique(phewas.res$id[i]))
    },
    error = function(e) {
      return(NA)
    }
  )
  
  # Merge with pQTL data.
  merge(select(pqtl.stats, -c(gene.snp.pair, ma_samples, ma_count, start_distance)),
        assoc.temp,
        by = c("chr", "position")
  ) %>%
    filter(!is.na(eaf)) %>%
    # Check the allele alignment.
    mutate(temp = variant_id) %>%
    separate(temp, into = c("chr", "pos", "ref.allele", "alt.allele"), sep = ":") %>%
    mutate(allele.check = ifelse(ea != ref.allele & ea != alt.allele, "inconsistent", "ok")) %>%
    filter(allele.check == "ok") %>%
    mutate(beta = ifelse(alt.allele == ea, beta, -1 * beta), MAF = ifelse(alt.allele != ea, 1 - eaf, eaf)) %>%
    mutate(top.qtl = rsID[which.max(abs(slope / slope_se))], top.gwas.snp = rsID[which.max(abs(beta / se))]) %>%
    arrange(position) %>%
    select(-c(chr, pos, allele.check)) -> pqtl.assoc
  
  if (nrow(pqtl.assoc) > 0) {
    # Store D1 D2 for each trait.
    D1 <- list(
      beta = pqtl.assoc$slope, varbeta = pqtl.assoc$slope_se^2,
      N = 191, sdY = 1, type = "quant",
      MAF = pqtl.assoc$af, snp = pqtl.assoc$rsid, position = pqtl.assoc$position
    )
    
    if (!("ncase" %in% names(trait.info))) {
      D2 <- list(
        beta = pqtl.assoc$beta, varbeta = pqtl.assoc$se^2,
        N = max(pqtl.assoc$n), type = "quant",
        MAF = pqtl.assoc$MAF, snp = pqtl.assoc$rsid, position = pqtl.assoc$position
      )
    } else {
      D2 <- list(
        beta = pqtl.assoc$beta, varbeta = pqtl.assoc$se^2,
        N = trait.info$sample_size, type = "cc",
        MAF = pqtl.assoc$MAF, snp = pqtl.assoc$rsid, position = pqtl.assoc$position,
        s = trait.info$ncase / (trait.info$ncontrol + trait.info$ncase)
      )
    }
    
    coloc.res <- coloc.signals(D1, D2, method = "single", p12 = 1e-5)

    coloc.res$summary$trait <- unique(assoc.temp$trait)
    coloc.res$summary$gene <- unique(assoc.temp$gene)
    coloc.res$summary$lead.snp <- unique(assoc.temp$lead.snp)
    coloc.res$summary$queried.snp <- unique(assoc.temp$queried.snp)
    coloc.res$summary$top.gwas.snp <- unique(pqtl.assoc$top.gwas.snp) # CHECK!
    
    if (coloc.res$summary$PP.H4.abf > 0) {
      return(list(pqtl.assoc = pqtl.assoc, coloc.res = coloc.res))
    } else {
      return(list(pqtl.assoc = pqtl.assoc, coloc.res = NULL))
    }
  } else {
    return(list(pqtl.assoc = pqtl.assoc, coloc.res = NULL))
  }
})

# Merge with the phewas results.
coloc.output <- lapply(coloc.res, function(x) {
  x$coloc.res$summary
})

coloc.output <- do.call(rbind, coloc.output)

# Extract also the SNP associations useful for regional plots.
pqtl.assoc <- lapply(coloc.res, function(x) {
  x$pqtl.assoc
})

pqtl.assoc <- do.call(rbind, pqtl.assoc)

# Remove columns.
to.select <- c("hit1.margz", "hit2.margz")

merged.df <- merge(phewas.res, select(coloc.output, -c(to.select)), by = c("trait", "queried.snp"))

# Export results.
out.prefix <- unique(pqtl.stats$context)

write.table(merged.df, paste("results/coloc.res", unique(merged.df$gene), out.prefix, unique(gsub("[:-]", ".", merged.df$region)), "txt", sep = "."), sep = "\t", row.names = F, quote = F)

write.table(pqtl.assoc, paste("results/pqtl.assoc", unique(pqtl.assoc$gene), out.prefix, "txt", sep = "."), sep = "\t", row.names = F, quote = F)

require("regional.plot.function.R")

setwd("")

# LBR.
regional.plot(
  assoc.df.pr1 = fread("results/pqtl.assoc.LBR.PR1.txt"), 
  assoc.df.pr2 = fread("results/pqtl.assoc.LBR.PR2.txt"),
  assoc.df.nr1 = fread("results/pqtl.assoc.LBR.NR1.txt"),
  assoc.df.nr2 = fread("results/pqtl.assoc.LBR.NR2.txt"),
  coloc.df = fread("results/coloc.res.LBR.PR2.1.224615737.226615737.txt"),
  phenotype = NULL,
  out.prefix = "PR2.stacked"
)

# MSRA.
regional.plot(
  assoc.df.pr1 = fread("results/pqtl.assoc.MSRA.PR1.txt"), 
  assoc.df.pr2 = fread("results/pqtl.assoc.MSRA.PR2.txt"),
  assoc.df.nr1 = fread("results/pqtl.assoc.MSRA.NR1.txt"),
  assoc.df.nr2 = fread("results/pqtl.assoc.MSRA.NR2.txt"),
  coloc.df = fread("results/coloc.res.MSRA.PR2.8.8911778.10911778.txt"),
  phenotype = NULL,
  out.prefix = "PR2.stacked"
)

# IL12RB1.
regional.plot(
  assoc.df.pr1 = fread("results/pqtl.assoc.IL12RB1.PR1.txt"), 
  assoc.df.pr2 = fread("results/pqtl.assoc.IL12RB1.PR2.txt"),
  assoc.df.nr1 = fread("results/pqtl.assoc.IL12RB1.NR1.txt"),
  assoc.df.nr2 = fread("results/pqtl.assoc.IL12RB1.NR2.txt"),
  coloc.df = fread("results/coloc.res.IL12RB1.NR1.19.17197813.19197813.txt"),
  phenotype = NULL,
  out.prefix = "NR1.stacked"
)

# PDLIM7.
regional.plot(
  assoc.df.pr1 = fread("results/pqtl.assoc.PDLIM7.PR1.txt"), 
  assoc.df.pr2 = fread("results/pqtl.assoc.PDLIM7.PR2.txt"),
  assoc.df.nr1 = fread("results/pqtl.assoc.PDLIM7.NR1.txt"),
  assoc.df.nr2 = fread("results/pqtl.assoc.PDLIM7.NR2.txt"),
  coloc.df = fread("results/coloc.res.PDLIM7.NR1.5.175924584.177924584.txt"),
  phenotype = NULL,
  out.prefix = "NR1.stacked"
)

# CASP3.
regional.plot(
  assoc.df.pr1 = fread("results/pqtl.assoc.CASP3.PR1.txt"), 
  assoc.df.pr2 = fread("results/pqtl.assoc.CASP3.PR2.txt"),
  assoc.df.nr1 = fread("results/pqtl.assoc.CASP3.NR1.txt"),
  assoc.df.nr2 = fread("results/pqtl.assoc.CASP3.NR2.txt"),
  coloc.df = fread("results/coloc.res.CASP3.NR1.4.184570629.186570629.txt"),
  phenotype = NULL,
  out.prefix = "NR1.stacked"
)

# FGF21.
regional.plot(
  assoc.df.pr1 = fread("results/pqtl.assoc.FGF21.PR1.txt"), 
  assoc.df.pr2 = fread("results/pqtl.assoc.FGF21.PR2.txt"),
  assoc.df.nr1 = fread("results/pqtl.assoc.FGF21.NR1.txt"),
  assoc.df.nr2 = fread("results/pqtl.assoc.FGF21.NR2.txt"),
  coloc.df = fread("results/coloc.res.FGF21.PR2.19.48258816.50258816.txt"),
  phenotype = "Eosinophil counts",
  out.prefix = "PR2"
)

regional.plot(
  assoc.df.pr1 = fread("results/pqtl.assoc.FGF21.PR1.txt"), 
  assoc.df.pr2 = fread("results/pqtl.assoc.FGF21.PR2.txt"),
  assoc.df.nr1 = fread("results/pqtl.assoc.FGF21.NR1.txt"),
  assoc.df.nr2 = fread("results/pqtl.assoc.FGF21.NR2.txt"),
  coloc.df = fread("results/coloc.res.FGF21.PR2.19.48258816.50258816.txt"),
  phenotype = "Plateletcrit",
  out.prefix = "PR2"
)
