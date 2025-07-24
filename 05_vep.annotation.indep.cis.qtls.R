# Load VEP annotations for all independent cis-pQTLs (n = ) across contexts.
vcf <- fread("results/indep.pqtls.vep.txt") # Results from online VEP query at 27/02/2025 - https://www.ensembl.org/Tools/VEP with default parameters.

# One variant can have >1 consequences which one should I keep?
# The consequence of the matched gene.
vcf %>% 
  mutate(IMPACT_score = ifelse(IMPACT == "HIGH", 4, 
                               ifelse(IMPACT=="MODERATE", 3, 
                                      ifelse(IMPACT == "LOW", 2,1) ))) -> vcf

# Retain the most significant IMPACT score.
vcf %>% 
  group_by(Location) %>% 
  dplyr::slice(which.max(IMPACT_score)) -> vcf

vcf %>% 
  group_by(Consequence) %>% 
  filter(n() > 2) %>% 
  select(Consequence) -> to.keep

# Which are the missense variants?
vcf %>% 
  filter(Consequence == "missense_variant")

vcf %>% 
  filter(Consequence %in% to.keep$Consequence) %>% 
  group_by(Consequence) %>% 
  summarise(fraction = round(n()/nrow(to.keep), 5)) -> plot.data

library(RColorBrewer)
getpal <- colorRampPalette(brewer.pal(9, "Set2"))
x <- getpal(20)

annotation.plot <- ggplot(data = plot.data,
                          aes(x = Consequence,
                              y = fraction,
                              fill = Consequence)) +
  theme_classic(base_size = 7) +
  geom_col(colour = "black") +
  scale_fill_manual(values = x) +
  labs(x = "Regulatory consequence", y = "Fraction of cis-pQTLs", fill = NULL) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/vep.indep.cis.pqtls.png", 
       annotation.plot, width = 12, height = 8.5, units = "cm", dpi = 300)