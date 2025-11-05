# -------------------------------------------------------------------------
# File: metaplot-shiny-intermediates
# Pierce et al 2025
# Purpose: Create metaplot intermediates for metaplots_shiny
# -------------------------------------------------------------------------

# Libraries ---------------------------------------------------------------
source("functions.R")

###metaplots

GBregion<-fread("../Output-files/Arabidopsis-PCSD-GeneBodyRegion.csv")
GBregion<-GBregion[totalintronNum!=0]
GBregion$IntronPct<-(GBregion$totalIntronLength/GBregion$geneLength)
epigmarks<-grep("H3|H2|H4|ATAC|IBM|SDG|Pol|JMJ", colnames(GBregion), value = TRUE)

featuresMetaplotsData <- function(data, group_by_column, grouping_column = "intronGp", marks, g = 10) {
  # Step 1: Create groups using cut2
  data[[grouping_column]] <- cut2(data[[group_by_column]], g = g)

  # Step 2: Summarize data (replace summarizeData with your actual function or summarize manually)
  summary <- summarizeData(data, marks, grouping_column)
}

# Define features and labels
features <- c("firstintronpos", "firstIntronLength", "totalintronNum", "totalIntronLength",
              "totalExonLength", "geneLength", "IntronPct")


epigmarks <- grep("H2|H4", colnames(GBregion), value = TRUE)
results1 <- lapply(seq_along(features), function(i) {
  df <- featuresMetaplotsData(GBregion, features[i], marks = epigmarks)
  df$geneFeature <- features[i]
  return(df)
})
combined1 <- rbindlist(results1)


epigmarks<-grep("ATAC|IBM|SDG|Pol|JMJ", colnames(GBregion), value = TRUE)
results2 <- lapply(seq_along(features), function(i) {
  df <- featuresMetaplotsData(GBregion, features[i], marks = epigmarks)
  df$geneFeature <- features[i]
  return(df)
})
combined2 <- rbindlist(results2)


epigmarks<-grep("H3K4|H3.1$|H3.3$", colnames(GBregion), value = TRUE)
results3 <- lapply(seq_along(features), function(i) {
  df <- featuresMetaplotsData(GBregion, features[i], marks = epigmarks)
  df$geneFeature <- features[i]
  return(df)
})
combined3 <- rbindlist(results3)


epigmarks<-grep("H3K14ac|H3K23ac|H3K27", colnames(GBregion), value = TRUE)
results4 <- lapply(seq_along(features), function(i) {
  df <- featuresMetaplotsData(GBregion, features[i], marks = epigmarks)
  df$geneFeature <- features[i]
  return(df)
})
combined4 <- rbindlist(results4)


epigmarks<-grep("H3K36|H3K56|H3K9", colnames(GBregion), value = TRUE)
results5 <- lapply(seq_along(features), function(i) {
  df <- featuresMetaplotsData(GBregion, features[i], marks = epigmarks)
  df$geneFeature <- features[i]
  return(df)
})
combined5 <- rbindlist(results5)


all<-rbind(combined1, combined2, combined3, combined4, combined5)
saveRDS(all, "metaplotsData.rds")

