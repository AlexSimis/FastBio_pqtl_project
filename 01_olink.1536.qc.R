## This notebook contains all the necessary steps to produce QCed NPX protein values from Olink 1536 explore.

library(tidyverse)
library(data.table)

# Loading data.

args <- commandArgs(trailingOnly = TRUE)
path = getwd()

if (length(args) < 2) {
  stop("Please provide two file paths: 1) NPX data file, 2) Metadata file.")
}

phenotype.data <- args[1]
meta.data <- args[2]

# Phenotype data.
raw.data <- as.data.frame(readr::read_delim(phenotype.data,
                                            delim = ";", 
                                            col_names =T, 
                                            show_col_types = F ))

# Meta data.
meta <- readr::read_delim(meta.data,
                          col_names=TRUE,
                          show_col_types = FALSE)

colnames(meta)[1] = 'SampleID' # Set proper column name.

# Remove control Samples.
qc.data <- filter(raw.data, !SampleID%in%c('SC1','SC2'))

# Duplicated Assays? - Keep the ones with the lowest missing frequencies.
qc.data %>% 
  distinct(OlinkID, .keep_all = T) %>% 
  group_by(Assay) %>% 
  filter(n()>1) %>% 
  slice(which.min(MissingFreq)) %>% 
  select(OlinkID) -> to.keep

qc.data %>% 
  distinct(OlinkID, .keep_all = T) %>% 
  group_by(Assay) %>% 
  filter(n()>1 & !OlinkID%in%to.keep$OlinkID) -> to.drop

qc.data <- filter(qc.data, !OlinkID%in%to.drop$OlinkID)

# Replace NPX values with NA when they have QC warning.
qc.data %>% 
  mutate(NPX = ifelse(QC_Warning =='WARN', NA, NPX)) -> qc.data

table(qc.data$QC_Warning, qc.data$Assay_Warning)

table(is.na(qc.data))

table(is.na(qc.data))[2]/(801*1472) # We are imputing about 1%.

# Conservative approach.
# Drop Assays with missing freq higher than 40%.
qc.data <- filter(qc.data, MissingFreq < 0.4)

length(unique(qc.data$OlinkID)) # Number of unique remained Assays.
length(unique(qc.data$Assay))

# Pivoting data to merge with metadata.
qc.data.wide <- tidyr::pivot_wider(dplyr::select(qc.data, SampleID, OlinkID, PlateID, NPX), names_from = 'OlinkID', values_from = 'NPX', values_fn = list)
qc.data.wide <- qc.data.wide %>% mutate_if(is.list, ~as.numeric(as.character(.))) # Convert the lists to numeric (dbl) 

length(select(qc.data.wide, starts_with("OID"))) # Sanity check.

table(is.na(qc.data.wide))

qc.data.wide.merged <- merge(meta, qc.data.wide, by = 'SampleID')

qc.data.wide.merged %>% 
  group_by(fasting, timePoint) %>% 
  summarise(n())

# Imputation of missing values (as produced by QC_Warning).
# Impute missing data via the mean per Assay for each group separately.
qc.data.wide.merged.imp <- qc.data.wide.merged

# Columns to apply the imputation. 
cols <- colnames(dplyr::select(qc.data.wide.merged.imp, starts_with('OID')))

# Context specific umputation.
qc.data.wide.merged.imp[which(qc.data.wide.merged.imp$fasting=='F' & qc.data.wide.merged.imp$timePoint=='T1'),  cols] <- apply(qc.data.wide.merged.imp[which(qc.data.wide.merged.imp$fasting=='F' & qc.data.wide.merged.imp$timePoint=='T1'),cols], 2, function(z){z = as.numeric(z); z[which(is.na(z))] = mean(z, na.rm = T); return(z)}) # for PV1
qc.data.wide.merged.imp[which(qc.data.wide.merged.imp$fasting=='F' & qc.data.wide.merged.imp$timePoint=='T2'),cols] <- apply(qc.data.wide.merged.imp[which(qc.data.wide.merged.imp$fasting=='F' & qc.data.wide.merged.imp$timePoint=='T2'),cols], 2, function(z){z = as.numeric(z); z[which(is.na(z))] = mean(z, na.rm = T); return(z)}) # for PV2
qc.data.wide.merged.imp[which(qc.data.wide.merged.imp$fasting=='NF' & qc.data.wide.merged.imp$timePoint=='T1'),cols] <- apply(qc.data.wide.merged.imp[which(qc.data.wide.merged.imp$fasting=='NF' & qc.data.wide.merged.imp$timePoint=='T1'),cols], 2, function(z){z = as.numeric(z); z[which(is.na(z))] = mean(z, na.rm = T); return(z)}) # for NV1
qc.data.wide.merged.imp[which(qc.data.wide.merged.imp$fasting=='NF' & qc.data.wide.merged.imp$timePoint=='T2'),cols] <- apply(qc.data.wide.merged.imp[which(qc.data.wide.merged.imp$fasting=='NF' & qc.data.wide.merged.imp$timePoint=='T2'),cols], 2, function(z){z = as.numeric(z); z[which(is.na(z))] = mean(z, na.rm = T); return(z)}) # for NV2

# Data to normalize.
transf.data <- qc.data.wide.merged.imp[,c("SampleID","fasting", "timePoint", cols)]
# Inverse normal transformation (INT) function.
inverseNormal <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}

# A function to apply the transformation to a dataframe.
transformation <- function(dataframe){
  normalised.df <- as.data.frame(apply(dataframe, 2, inverseNormal))
  normalised.df
}

# At this point the transformation could take place per group of interest or globally.
transf.data[,cols] <- transformation(transf.data[,cols])

# Export the QCed data.
out.prefix = "QCed&transformed.data.txt"
path = "testing/"
caroline::write.delim(transf.data,paste0(path, out.prefix), quote = F)

# Sanity check.
#temp <- fread(paste0(path, "QCed&transformed.data.txt"))

