# -------------------------------------------------------------------------
# File: supplementary-tables
# Pierce et al 2025
# Purpose: Create supplementary tables
# -------------------------------------------------------------------------

# Libraries ---------------------------------------------------------------
source("functions.R")

# File Correlations -------------------------------------------------------

# List all CSV files in a directory
file_list <- list.files(path = "../Output-files/file-correlations/", pattern = "\\.csv$", full.names = TRUE)

# Read them all into a list of data.tables
data_list <- lapply(file_list, fread)

# Step 2: Read and reshape each file
long_corr_list <- lapply(file_list, function(file) {
  # Read the correlation table
  dt <- fread(file)

  # Melt the table to long format (file1 = row, file2 = column)
  long_dt <- melt(dt, id.vars = "Variable", variable.name = "File 2", value.name = "corr")

  # Rename columns
  setnames(long_dt, "Variable", "File 1")

  # Optionally add a source column
  #long_dt[, source_file := basename(file)]

  return(long_dt)
})

# Step 3: Combine all into one data.table
all_long_corr <- rbindlist(long_corr_list)

# Convert factors to characters
all_long_corr[, `File 1` := as.character(`File 1`)]
all_long_corr[, `File 2` := as.character(`File 2`)]

# Generate canonical pairs
all_long_corr[, `:=`(
  file_a = pmin(`File 1`, `File 2`),
  file_b = pmax(`File 1`, `File 2`)
)]

# Remove duplicates
unique_corrs <- all_long_corr[`File 1` == file_a & `File 2` == file_b, .(`File 1`, `File 2`, corr)]
unique_corrs$Mark<-gsub("-.*", "", unique_corrs$`File 1`)
unique_corrs <- unique_corrs[, .(Mark, `File 1`, `File 2`, corr)]
write.csv(unique_corrs, "../Output-files/File-Correlations.csv", row.names = FALSE)

# Spearman Rho ------------------------------------------------------------

# List all CSV files in a directory
spearmanRho_list <- list.files(path = "../Output-files/CorrelationStatistics/", pattern = "\\Rho.csv$", full.names = TRUE)

# Extract just the base file names (without path)
file_names <- basename(spearmanRho_list)

# Read and name the list
Rho_list <- setNames(lapply(spearmanRho_list, fread), file_names)

# Add the file name as a new column in each data.table
Rho_list_named <- Map(function(dt, name) {
  dt[, `Gene Region` := name]
  return(dt)
}, Rho_list, names(Rho_list))

# Combine into a single data.table
Rho_all <- rbindlist(Rho_list_named)
Rho_all$`Gene Region`<-gsub("Rho.csv", "", Rho_all$`Gene Region`)
Rho_all$`Gene Region`<-gsub("AllGene", "GeneSequence", Rho_all$`Gene Region`)
Rho_all$`Gene Region`<-gsub("CDS", "ExonSequence", Rho_all$`Gene Region`)
Rho_all$`Gene Region`<-gsub("NCS", "IntronSequence", Rho_all$`Gene Region`)

setcolorder(Rho_all, c("Gene Region", setdiff(names(Rho_all), "Gene Region")))
write.csv(Rho_all, "../Output-files/AllRho.csv", row.names = FALSE)

# Spearman P-Value --------------------------------------------------------

# List all CSV files in a directory
spearmanPVal_list <- list.files(path = "../Output-files/CorrelationStatistics/", pattern = "\\PVal.csv$", full.names = TRUE)

# Extract just the base file names (without path)
file_names <- basename(spearmanPVal_list)

# Read and name the list
PVal_list <- setNames(lapply(spearmanPVal_list, fread), file_names)

# Add the file name as a new column in each data.table
PVal_list_named <- Map(function(dt, name) {
  dt[, `Gene Region` := name]
  return(dt)
}, PVal_list, names(PVal_list))

# Combine into a single data.table
PVal_all <- rbindlist(PVal_list_named)
PVal_all$`Gene Region`<-gsub("PVal.csv", "", PVal_all$`Gene Region`)
PVal_all$`Gene Region`<-gsub("AllGene", "GeneSequence", PVal_all$`Gene Region`)
PVal_all$`Gene Region`<-gsub("CDS", "ExonSequence", PVal_all$`Gene Region`)
PVal_all$`Gene Region`<-gsub("NCS", "IntronSequence", PVal_all$`Gene Region`)

setcolorder(PVal_all, c("Gene Region", setdiff(names(PVal_all), "Gene Region")))
write.csv(PVal_all, "../Output-files/AllPValue.csv", row.names = FALSE)

# Model Validation --------------------------------------------------------
### Metrics
AllGeneMetrics<-fread("../Output-files/Arabidopsis-PCSD-CVMetrics.csv")
AllGeneMetrics$`Gene Region`<-"GeneSequence"
ExonMetrics<-fread("../Output-files/Arabidopsis-PCSD-CDSCVMetrics.csv")
ExonMetrics$`Gene Region`<-"ExonSequence"
IntronMetrics<-fread("../Output-files/Arabidopsis-PCSD-NonCodingCVMetrics.csv")
IntronMetrics$`Gene Region`<-"IntronSequence"

AllMetrics<-rbind(AllGeneMetrics, ExonMetrics, IntronMetrics)
setcolorder(AllMetrics, c("Gene Region", "mark", "RMSE", "MAE", "R2", "COR"))
write.csv(AllMetrics, "../Output-files/AllMetrics.csv", row.names = FALSE)

### Features
AllGeneFeatures<-fread("../Output-files/Arabidopsis-PCSD-CVFeatures.csv")
AllGeneFeatures$`Gene Region`<-"GeneSequence"
ExonFeatures<-fread("../Output-files/Arabidopsis-PCSD-CDSCVFeatures.csv")
ExonFeatures$`Gene Region`<-"ExonSequence"
IntronFeatures<-fread("../Output-files/Arabidopsis-PCSD-NonCodingCVFeatures.csv")
IntronFeatures$`Gene Region`<-"IntronSequence"

AllFeatures<-rbind(AllGeneFeatures, ExonFeatures, IntronFeatures)
setcolorder(AllFeatures, c("Gene Region", "mark", "rn", "%IncMSE", "IncNodePurity", "CV"))
write.csv(AllFeatures, "../Output-files/AllFeatures.csv", row.names = FALSE)


# Excel -------------------------------------------------------------------
FileCorrelations<-fread("../Output-files/File-Correlations.csv")
AllRho<-fread("../Output-files/AllRho.csv")
AllPValue<-fread("../Output-files/AllPValue.csv")
Metrics<-fread("../Output-files/AllMetrics.csv")
Features<-fread("../Output-files/AllFeatures.csv")

data_list <- list(
  "Supplementary Table S1" = FileCorrelations, # ChIP-seq and ATAC-seq file correlations
  "Supplementary Table S2" = AllRho, # Pairwise Spearman Correlation Rho between intron architecture and chromatin features
  "Supplementary Table S3" = AllPValue, # Statistical Significance (p-values) of Spearman Correlations between intron architecture and chromatin features
  "Supplementary Table S4" = Metrics, # Random Forest Model Performance Metrics
  "Supplementary Table S5" = Features # Random Forest Model Feature Importance
)

# Write to Excel file
write_xlsx(data_list, "../Output-files/SupplementaryTables.xlsx")
