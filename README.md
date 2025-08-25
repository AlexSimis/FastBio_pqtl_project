# cis-pQTL analysis in the FastBio study
+ Link to preprint in BioRxiv: https://www.biorxiv.org/content/10.1101/2025.07.21.665884v1
+ Link to summary statistics produced in this study: pending

# Contents
## Files:
+ ***01_olink.1536.qc.R*** script containing all the necessary steps for processing and normalizing the raw proteomic data (NPX values) from Olink 1536 Explore
+ ***02_genotype.pc.estimation.sh*** shell script used to estimate the genotype-derived principal components using plinkv2. These were used as covariates in the main cis-pQTL analysis.
+ ***03_covariate.selection.R*** script used to select the covariates to adjust for in the main state-state cis-pQTL analysis. These include both fixed covariates, such as age, sex, BMI and medication_use, as well as latent factors derived from PC analysis of protein and genotype data.
+ ***04_tensor.cis.qtl.mapping.sh*** shell script used to perform the following tasks using TensorQTL: (a) Perform cis-pQTL mapping with permutations and default parameters, (b) full summary statistics per protein tested, (c) independent cis-pQTL mapping using default parameters.
+ ***05_vep.annotation.indep.cis.qtls.R*** script to explore the annotations of all independent cis-pQTLs across contexts. The annotations were retrieved using the VEP online tool (https://www.ensembl.org/Tools/VEP)
+ ***06_top.cis.qtl.replication.R*** script used to explore the replication of top cis-pQTLs against the cis-summary statistics of the latest UKBB-PPP publication. Data were freely accessed using the following link https://www.synapse.org/Synapse:syn51365303.
+ ***07_mashr.cis.qtl.mapping.R*** script used to jointly analyze the results of steady-state cis-pQTLs from all contexts and decipher which SNP-protein pairs exhibit differential effects between contexts.
+ ***08_colocalization.analysis.R*** script used to (a) perform phewas queries of differential cis-pQTLs against OpenGWAS catalog and (b) based on the associated traits run colocalization per context.
+ ***regional.plot.function.R*** this function was used to produce stacked regional plots based on the colocalization results across contexts.
