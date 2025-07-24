# This is a function to produce stacked regional plots of colocalization results per context, using "locuszoomr"

regional.plot <- function(assoc.df.pr1, assoc.df.pr2, assoc.df.nr1, assoc.df.nr2, coloc.df, phenotype, out.prefix) {
  require(data.table)
  require(tidyverse)
  if(is.null(phenotype)){
    # Find which trait has the highest PPH4.
    coloc.df %>%
      dplyr::slice(which.max(PP.H4.abf)) -> coloc.info
  }  else {
    coloc.df %>% 
      dplyr::filter(trait == phenotype) -> coloc.info
  }
  
  # Find the context of sig. coloc.
  bind_rows(assoc.df.pr1, assoc.df.pr2, assoc.df.nr1, assoc.df.nr2) %>% 
    dplyr::filter(trait == coloc.info$trait) %>% 
    dplyr::slice(which.min(pval_nominal))  -> to.select
  
  # Filter the assoc file based on the phenotype-trait.
  rbind(assoc.df.pr1, assoc.df.pr2, assoc.df.nr1, assoc.df.nr2)%>%
    dplyr::filter(trait == coloc.info$trait & context == to.select$context) %>%
    mutate(temp = hg19) %>% # Use the hg19 build.
    separate(temp, into = c("chr", "pos"), sep = ":") -> temp
  
  
  # Store the pQTL stats as locus.data
  assoc.df.ref <- function(data){
    data %>% 
      mutate(temp = hg19) %>% # Use the hg19 build.
      separate(temp, into = c("chr", "pos"), sep = ":") %>% 
      dplyr::filter(hg19 %in% temp$hg19) %>% # To match the locus.data. %>% 
      mutate(
        chrom = chr, 
        pos = as.numeric(pos),
        rsid = rsID, 
        other_allele = ref.allele,
        effect_allele =  alt.allele,
        p = pval_nominal, 
        beta = slope, 
        se = slope_se) %>%
      dplyr::select(chrom, pos, rsid, other_allele, effect_allele, p, beta, se, phenotype_id, trait, context) %>%
      arrange(pos) 
    
  }
  
  locus.data.pr1 <-  assoc.df.ref(assoc.df.pr1)
  locus.data.pr2 <-  assoc.df.ref(assoc.df.pr2)
  locus.data.nr1 <-  assoc.df.ref(assoc.df.nr1)
  locus.data.nr2 <-  assoc.df.ref(assoc.df.nr2)
  
  # Store the trait stats as locus.data.gwas
  temp %>%
    mutate(
      chrom = chr, 
      pos = as.numeric(pos), 
      rsid = rsID,
      other_allele = ref.allele,
      effect_allele = alt.allele,
      p = p, 
      beta = beta, 
      se = se
    ) %>%
    dplyr::select(chrom, pos, rsid, other_allele, effect_allele, p, beta, se) %>% 
    arrange(pos) -> locus.data.gwas
  
  # Define locus data for locuszoomr based on the TSS of the candidate protein.
  # Window 2 Mb.
  chr <- as.numeric(unique(temp$chr))
  tss <- as.numeric(coloc.info$transcription_start_site)
  xrange <- c(tss - 1e6, tss + 1e6)
  phenotype_id <- unique(temp$phenotype_id)
  lead_snp <- to.select$rsID
  
  temp.prefix <- unique(temp$trait)
  # temp.prefix <- "Obesity"
  label_title <-  temp.prefix
  
  library(locuszoomr)
  require(EnsDb.Hsapiens.v75)
  
  # Create loci with LD.
  create_locus <- function(data) {
    locus(data = data, seqname = chr, xrange = xrange,
          ens_db = "EnsDb.Hsapiens.v75", LD = "r2") %>%
      link_LD(token = "335da46cccf2")
  }
  
  loc.pr1 <- create_locus(locus.data.pr1)
  loc.pr2 <- create_locus(locus.data.pr2)
  loc.nr1 <- create_locus(locus.data.nr1)
  loc.nr2 <- create_locus(locus.data.nr2)
  loc.gwas <- create_locus(locus.data.gwas)
  
  # Scatter plot function.
  make_scatter_plot <- function(loc, title, legend_pos) {
    gg_scatter(
      loc,
      index_snp = lead_snp,
      labels = NULL,
      legend_pos = NULL,
      beta = "beta",
      pcutoff = 1e-5,
      xlab = NULL,
      xticks = FALSE,
      size = 1,
      # shape_values = c(21, 24, 25)
    ) +
      geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "grey") +
      scale_y_continuous(limits = c(0, max(-log10(to.select$pval_nominal)))) +
      labs(title = title) +
      theme_classic(base_size = 7) +
      theme(legend.position="none")
    
  }
  
  # Build regional plots per context.
  pr1 <- make_scatter_plot(loc.pr1, "PR at T1", legend_pos = NULL)
  pr2 <- make_scatter_plot(loc.pr2, "PR at T2", legend_pos = NULL)
  nr1 <- make_scatter_plot(loc.nr1, "NR at T1", legend_pos = NULL)
  nr2 <- make_scatter_plot(loc.nr2, "NR at T2", legend_pos = NULL)
  
  # Make regional plot for GWAS.
  gwas <- gg_scatter(
    loc.gwas,
    # labels = c("index", lead_snp),
    labels = NULL,
    legend_pos = NULL,
    beta = "beta",
    pcutoff = 1e-5,
    xlab = NULL,
    xticks = FALSE,
    size = 1
    # shape_values = c(21, 24, 25)
  ) +
    labs(title = label_title) +
    theme_classic(base_size = 7) +
    theme(legend.position="none")
  
  
  # Add gene annotations.
  genes <- gg_genetracks(loc.gwas, 
                         highlight = phenotype_id, 
                         # filter_gene_name = phenotype_id,
                         maxrows = 4,
                         filter_gene_biotype = 'protein_coding',
                         cex.text = 0.4) +
    theme_classic(base_size = 7) +
    theme(legend.position="none")
  
  require(patchwork)
  comb.plot <- wrap_plots(pr1, pr2, nr1, nr2, gwas, genes, ncol = 1)
  
  require(ggpubr)
  
  # comb.plot <- comb.plot & common_theme
  comb.plot <- comb.plot & theme(text = element_text(size = 7) )
  
  assign(paste0("regional.plot.", out.prefix), comb.plot, envir = .GlobalEnv)
  
  
  # dev.new(width = 120, height = 180, unit = "mm", noRStudioGD=TRUE)
  # ggsave(paste0("results/plots/", paste(unique(temp$phenotype_id),
  #                                       temp.prefix,
  #                                       out.prefix,
  #                                       sep = "."), ".png"), comb.plot, width = 100, height = 180, unit = "mm", dpi = 300, scaling = 1)
  
  lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
  
}

