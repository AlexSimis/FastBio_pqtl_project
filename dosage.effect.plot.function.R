# Dosage plot function.
dosage.plot <- function(snp.matrix, variant.id, pheno1, pheno2, p1, p2, ylab, trait, out.prefix, colors, plot.title){
  # Extract the dosage effect of the selected SNP.
  dosage <- snps[, which(colnames(snps) == variant.id)]
  dosage <- data.frame(sampleID = names(dosage), dosage = as.matrix(dosage)[ ,1] )
  # dosage |> mutate(dosage = 2 - as.numeric(dosage)) -> dosage # This is used selectively e.g. FGF21 case.
  
  # Properly format the sampleID.
  sampleID <- str_split(dosage$sampleID, "_", simplify = T)
  sampleID <- paste(sampleID[, 1], sampleID[, 2], sep = "_")
  dosage$sampleID <- sampleID

  # Merge the dosage effect with the phenotype data.
  temp <- rbind(pheno1, pheno2)
  geno.pheno <- merge(dosage, select(temp, IID, trait), by.x = "sampleID", by.y = "IID")
  geno.pheno |> mutate(timepoint = gsub(".*_", "", sampleID)) |>
    mutate(temp = gsub("_.*", "", sampleID))  |>
    group_by(temp) |> # This is to ensure that we have the same sampleID for both time points.
    filter(n() > 1) -> geno.pheno
  
  print(table(geno.pheno$dosage))
  
  library(ggplot2)
  require(ggpubr)
  
  # Perform linear regressions and extract p-values for T1 and T2 time points
  lm1 <- summary(lm(as.formula(paste(trait, "~ dosage")), data = filter(geno.pheno, timepoint == "T1")))
  lm1.res <- lm1$coefficients[2,]
  lm1.pval <- lm1$coefficients[2, 4]
  
  lm2 <- summary(lm(as.formula(paste(trait, "~ dosage")), data = filter(geno.pheno, timepoint == "T2")))
  lm2.res <-  lm2$coefficients[2,]
  lm2.pval <- lm2$coefficients[2, 4]
  
  # Save lm. results to the global environment.
  lm.res <- rbind(lm1.res, lm2.res)
  assign("lm.res", lm.res, envir = .GlobalEnv)
  
  print(lm1); print(lm2)
  print(lm1.pval); print(lm2.pval)
  
  # get_asterisks <- function(p) {
  #   if (p < 0.05) return("bold")
  #   else return("")  # Not significant
  # }

  # T1 label
  if (is.null(p1)) {
    pt1 <- paste0(
      "p = ", formatC(lm1.pval, format = "e", digits = 2)
    )
    p1_font <- NULL
  } else {
    pt1 <- paste0(
      "p = ", formatC(lm1.pval, format = "e", digits = 2)
    )
    p1_font <- ifelse(p1 < 0.05, "bold", "plain")
  }
  
  # T2 label
  if (is.null(p2)) {
    pt2 <- paste0(
      "p = ", formatC(lm2.pval, format = "e", digits = 2)
      )
    p2_font <- NULL
  } else {
    pt2 <- paste0(
      "p = ", formatC(lm2.pval, format = "e", digits = 2)
    )
    p2_font <- ifelse(p2 < 0.05, "bold", "plain")
  }
  
  # This is hard coded for the specific case of validation of LBR diff. pqtl with cholesterol classes were we have used adjustment correction.
  # p1 <- "ns"
  # p2 <- "*"
  
  # Create the plot and add p-value annotations
  plot <- ggplot(data = geno.pheno, aes(x = as.factor(dosage), y = .data[[trait]], fill = timepoint)) +
    geom_violin(size = 0.5, trim = FALSE, scale = "width") +
    geom_jitter(size = 0.7, alpha = 0.4, width = 0.1) +
    geom_boxplot(width = 0.3, fill = "white", outlier.shape = NA, alpha = 0.5) +
    #scale_fill_brewer(palette = "Pastel1") +  # Use a color palette for the fill
    theme_classic(base_size = 7) +
    scale_fill_manual(values = colors, guide = "none") +
    labs(x = paste(variant.id, sep = " "), y = ylab, fill = "Timepoint") +  # Add legend for timepoint
    # Annotate p-value for T1
    geom_text(data = subset(geno.pheno, timepoint == "T1"), 
              aes(x = 1, y = 4.5, # max(abs(geno.pheno[[trait]]), na.rm = TRUE) * 1.6
                  label = pt1,
                  fontface = p1_font),
              size = 2.2, hjust = 0, vjust = 0, check_overlap = TRUE) +
    # Annotate p-value for T2
    geom_text(data = subset(geno.pheno, timepoint == "T2"), 
              aes(x = 1, y = 4.5, #4.5
                  label = pt2,
                  fontface = p2_font),
              size = 2.2,  hjust = 0, vjust = 0, check_overlap = TRUE) +
    labs(title = plot.title)  +
    theme(plot.title = element_text(size = 7, hjust = 0.5),
          axis.ticks = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid = element_line(color = "white"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.background = element_rect(fill = "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_x_discrete(labels = c(paste0("0\nN=",sum(geno.pheno$dosage == 0)/2), 
                                paste0("1\nN=",sum(geno.pheno$dosage == 1)/2), 
                                paste0("2\nN=",sum(geno.pheno$dosage == 2)/2))) +
    facet_wrap(~timepoint, scales = "fixed") +
    ylim(-3.5, 5.5)

  assign(out.prefix, plot, envir =.GlobalEnv)
  #ggsave(paste0(path, "dosage.plot.", variant.id, ".", out.prefix, ".png"), plot, width = 6.7, height = 6.1, units = "cm", dpi = 300)
}


