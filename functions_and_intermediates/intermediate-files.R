# -------------------------------------------------------------------------
# File: intermediate-files
# Pierce et al 2025
# Purpose: Create intermediate files for data analysis
# -------------------------------------------------------------------------

# Libraries ---------------------------------------------------------------
source("functions.R")

# Gene, Exon, Intron Coordinates & Gene Features table --------------------
# Data from TAIR10
ArabidopsisFeatures<-intronFeatures("../Input-files/TAIR10_GFF3_genes.gff")
write.csv(ArabidopsisFeatures$allGenes, "../Output-files/ArabidopsisGeneCoordinates.csv", row.names = FALSE)
write.csv(ArabidopsisFeatures$allExons, "../Output-files/ArabidopsisExonCoordinates.csv", row.names = FALSE)
write.csv(ArabidopsisFeatures$allIntrons, "../Output-files/ArabidopsisIntronCoordinates.csv", row.names = FALSE)
write.csv(ArabidopsisFeatures$IntronFeatures, "../Output-files/ArabidopsisGeneFeatures.csv", row.names = FALSE)

# 10 bp file --------------------------------------------------------------
# Data from PCSD
ArabidopsisScores<-epigenomeOverlaps("../BW-files/Arabidopsis-PCSD/")
write.csv(ArabidopsisScores, "../Output-files/Arabidopsis-PCSD-10.csv", row.names = FALSE)

# correlations between files for the same chromatin feature --------------
get_prefix <- function(colname) {
  strsplit(colname, "-")[[1]][1]
}
column_groups <- split(names(ArabidopsisScores[,-c("V1", "chr", "start", "stop", "window_id")]),
                       sapply(names(ArabidopsisScores[,-c("V1", "chr", "start", "stop", "window_id")]), get_prefix))

generate_correlation("ATAC")
generate_correlation("H2A.W")
generate_correlation("H2A.X")
generate_correlation("H2A.Z")
generate_correlation("H3.1")
generate_correlation("H3.3")
generate_correlation("H3K14ac")
generate_correlation("H3K23ac")
generate_correlation("H3K27ac")
generate_correlation("H3K27me1")
generate_correlation("H3K27me3")
generate_correlation("H3K36ac")
generate_correlation("H3K36me3")
generate_correlation("H3K4me1")
generate_correlation("H3K4me2")
generate_correlation("H3K4me3")
generate_correlation("H3K56ac")
generate_correlation("H3K9ac")
generate_correlation("H3K9me1")
generate_correlation("H3K9me2")
generate_correlation("H4K16ac")
generate_correlation("IBM1")
generate_correlation("JMJ14")
generate_correlation("Pol2")
generate_correlation("PolIV")
generate_correlation("SDG8")

# scale to mean of 0 -----------------------------------------------------
columns_to_scale <- colnames(ArabidopsisScores[, -c(1:4)])
ArabidopsisScores[, columns_to_scale] <- lapply(ArabidopsisScores[, ..columns_to_scale], scale)
write.csv(ArabidopsisScores, "../Output-files/Arabidopsis-PCSD-scaled-10.csv", row.names = FALSE)

# normalized sum ---------------------------------------------------------
normalized<-normalizedSum(ArabidopsisScores)
write.csv(normalized, "../Output-files/Arabidopsis-PCSD-norm-10.csv", row.names = FALSE)

# Gene Depth -------------------------------------------------------------
GeneScores<-geneScores(ArabidopsisFeatures$allGenes, normalized, "NormalizedSum")
GeneScores<-isRepresentative(GeneScores, "../Input-files/TAIR10_representative_gene_models.txt")
write.csv(GeneScores, "../Output-files/Arabidopsis-PCSD-GeneDepth.csv", row.names = FALSE)

# Merge Gene Depth with Intron Features ----------------------------------
GeneScores<-merge(ArabidopsisFeatures$IntronFeatures, GeneScores, by=c("gene", "chr", "strand"))
write.csv(GeneScores,"../Output-files/Arabidopsis-PCSD-GeneFeatures.csv", row.names = FALSE)

# All Intron & Exon Depth ---------------------------------------------------
IntronExonScores<-IntronExonScores(ArabidopsisFeatures$allExons, ArabidopsisFeatures$allIntrons, normalized, "NormalizedSum")
IntronExonScores<-isRepresentative(IntronExonScores, "../Input-files/TAIR10_representative_gene_models.txt")
write.csv(IntronExonScores, "../Output-files/Arabidopsis-PCSD-IntronExonDepth.csv", row.names = FALSE)

# CDS & Non-Coding Depth -----------------------------------

CDSNCTable <- IntronExonScores %>%
  group_by(gene, chr, strand,Feature) %>%
  summarise(
    FeatureLength = sum(Length),
    across(starts_with("NormalizedSum"), ~ sum(.x * Length, na.rm = TRUE) / sum(Length, na.rm = TRUE), .names = "weighted_mean_{.col}")
  )

write.csv(CDSNCTable, "../Output-files/Arabidopsis-PCSD-CDSandNCDepth.csv", row.names = FALSE)

# Expression Data -------------------------------------------------------
# Download data from https://www.nature.com/articles/s41586-020-2094-2
expn<-fread("../Input-files/expressionfile") #change this to local file
expn<-select(expn, gene, matches("^(adult|callus|cell|dry|embryo|flower|seed|silique)"))
expn[, (names(expn)) := lapply(.SD, function(x) replace(x, is.na(x), 0))] # replace NAs with 0
expn[, avgExpression := sum(.SD, na.rm = TRUE)/54, .SDcols = 2:55, by = 1:nrow(expn)]
expn[, tissue_breadth := rowSums(.SD != 0), .SDcols = 2:55, by = 1:nrow(expn)]
expn[, cv := apply(.SD, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)), .SDcols = 2:55, by = 1:nrow(expn)]
write.csv(expn, "../Output-files/Arabidopsis-tissueExpression.csv", row.names = FALSE)

# Methylation Data ------------------------------------------------------
# Follow instructions from: https://www.nature.com/articles/s41586-021-04269-6
fasta_index <- FaFile("../FASTA/TAIR10_chr_all.fas")
open(fasta_index)

methylation <- fread("../Input-files/methylationfile") #change this to local file
colnames(methylation)<-c("chr", "start", "stop", "strand", "V5", "V6", "Context", "Nts")
methylation$Value<-1
setkey(methylation, chr, strand, start, stop)

ArabidopsisFeatures$allGenes$chr<-paste0("Chr", ArabidopsisFeatures$allGenes$chr)
ArabidopsisFeatures$allIntrons$chr<-paste0("Chr", ArabidopsisFeatures$allIntrons$chr)
ArabidopsisFeatures$allExons$chr<-paste0("Chr", ArabidopsisFeatures$allExons$chr)

# genes
wide<-extractCytme(fasta_index, ArabidopsisFeatures$allGenes, methylation)
write.csv(wide, "../Output-files/Arabidopsis-methylation.csv", row.names = FALSE)

#introns
intronCytme<-extractCytmeIntron(fasta_index, ArabidopsisFeatures$allIntrons, methylation)
write.csv(intronCytme, "../Output-files/Arabidopsis-intronmethylation.csv", row.names = FALSE)

# Calculate gene-level methylation across all introns
gene_level_methylation <- intronCytme %>%
  group_by(gene) %>%
  summarise(
    TotalCytCount = sum(CytCount, na.rm = TRUE),
    TotalmeCount_CG = sum(CG * CytCount, na.rm = TRUE),  # Methylation count for CG
    TotalmeCount_CHG = sum(CHG * CytCount, na.rm = TRUE), # Methylation count for CHG
    TotalmeCount_CHH = sum(CHH * CytCount, na.rm = TRUE), # Methylation count for CHH
    # Methylation ratio across all introns
    GeneMethylationRatio_CG = TotalmeCount_CG / TotalCytCount,
    GeneMethylationRatio_CHG = TotalmeCount_CHG / TotalCytCount,
    GeneMethylationRatio_CHH = TotalmeCount_CHH / TotalCytCount
  )
write.csv(gene_level_methylation, "../Output-files/Arabidopsis-noncodingmethylation.csv", row.names = FALSE)

#exons
exonCyt<-extractCytmeExon(fasta_index, ArabidopsisFeatures$allExons, methylation)
write.csv(exonCyt, "../Output-files/Arabidopsis-exonmethylation.csv", row.names = FALSE)

gene_levelexon_methylation <- exonCyt %>%
  group_by(gene) %>%
  summarise(
    TotalCytCount = sum(CytCount, na.rm = TRUE),
    TotalmeCount_CG = sum(CG * CytCount, na.rm = TRUE),  # Methylation count for CG
    TotalmeCount_CHG = sum(CHG * CytCount, na.rm = TRUE), # Methylation count for CHG
    TotalmeCount_CHH = sum(CHH * CytCount, na.rm = TRUE), # Methylation count for CHH
    # Methylation ratio across all exons
    GeneMethylationRatio_CG = TotalmeCount_CG / TotalCytCount,
    GeneMethylationRatio_CHG = TotalmeCount_CHG / TotalCytCount,
    GeneMethylationRatio_CHH = TotalmeCount_CHH / TotalCytCount
  )
write.csv(gene_levelexon_methylation, "../Output-files/Arabidopsis-CDSmethylation.csv", row.names = FALSE)

# Master Data Sets ----------------------------------------------

# All Gene
expn<-select(expn, gene, avgExpression, tissue_breadth, cv)
GeneScores$gene<-gsub("\\..*","", GeneScores$gene)

MasterData<-merge(GeneScores, expn, by="gene")

wide<-isRepresentative(wide,  "../Input-files/TAIR10_representative_gene_models.txt")
wide$gene<-gsub("\\..*","", wide$gene)
MasterData<-merge(MasterData, wide, by="gene")
MasterData<-MasterData[,ExonNum:=NULL]
MasterData$IntronPct<-(MasterData$totalIntronLength/MasterData$geneLength)

write.csv(MasterData, "../Output-files/Arabidopsis-PCSD-MasterData.csv", row.names = FALSE)

# CDS Only
CDSNCTable$gene<-gsub("\\..*","", CDSNCTable$gene)
CDSNCTable<-merge(CDSNCTable, expn, by="gene")
IntronPct<-MasterData %>% select(gene, IntronPct, firstIntronLength, firstintronpos,
                                 geneLength, totalExonLength, totalIntronLength, totalintronNum)

CDSNCTable<-as.data.table(CDSNCTable)
ExonScores<-CDSNCTable[Feature=="Exon"]
exonMethylation<-fread("../Output-files/Arabidopsis-CDSmethylation.csv")
exonMethylation<-isRepresentative(exonMethylation, "../Input-files/TAIR10_representative_gene_models.txt")
exonMethylation$gene<-gsub("\\..*","", exonMethylation$gene)
ExonScores<-merge(ExonScores, exonMethylation, by="gene")
ExonScores<-merge(ExonScores, IntronPct, by="gene")
write.csv(ExonScores, "../Output-files/Arabidopsis-PCSD-ExonMasterData.csv", row.names = FALSE)

# Non-coding Only
IntronScores<-CDSNCTable[Feature=="Intron"]
intronMethylation<-fread("../Output-files/Arabidopsis-noncodingmethylation.csv")
intronMethylation<-isRepresentative(intronMethylation, "../Input-files/TAIR10_representative_gene_models.txt")
intronMethylation$gene<-gsub("\\..*","", intronMethylation$gene)
IntronScores<-merge(IntronScores, intronMethylation, by="gene")
IntronScores<-merge(IntronScores, IntronPct, by="gene")
write.csv(IntronScores, "../Output-files/Arabidopsis-PCSD-IntronMasterData.csv", row.names = FALSE)

# Genomewide correlations --------------------------------------------

### All Gene
setnames(MasterData, old = names(MasterData), new = sub("^NormalizedSum_", "", names(MasterData)))
SubsetMD<-select(MasterData, starts_with(c("H2","H3", "ATAC", "H4", "IBM", "JMJ",
                                           "Pol", "SDG", "avg", "tissue", "cv",
                                           "CG", "CHH","CHG", "first", "total","Intron")), geneLength)

epigmarks <- grep("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC|avg|tissue|cv|IBM|JMJ|Pol|SDG", colnames(SubsetMD), value = TRUE)

# for all genes
corAllGenes<-correlationDataTable(SubsetMD, epigmarks, c("totalintronNum", "geneLength", "totalExonLength"))

# for intron-containing genes
corIntronGenes<-correlationDataTable(SubsetMD[totalintronNum!=0], epigmarks,
                                     c("firstIntronLength", "firstintronpos",
                                       "totalIntronLength", "IntronPct"))

correlations<-rbind(corAllGenes, corIntronGenes)
write.csv(correlations, "../Output-files/Arabidopsis-GWSpearmanCor.csv", row.names = FALSE)

rhoAndPValAllGenes<-correlationDataTableWithPValues(SubsetMD, epigmarks, c("totalintronNum", "geneLength", "totalExonLength"))
RhoTable<- rhoAndPValAllGenes$RhoTable %>% pivot_wider(names_from = Var1, values_from = Rho)
PValueTable<- rhoAndPValAllGenes$PValueTable %>% pivot_wider(names_from = Var1, values_from = PValue)

rhoAndPValIntronGenes<-correlationDataTableWithPValues(SubsetMD[totalintronNum!=0], epigmarks,
                                                       c("firstIntronLength", "firstintronpos",
                                                         "totalIntronLength", "IntronPct"))
RhoIntronTable<-rhoAndPValIntronGenes$RhoTable %>% pivot_wider(names_from = Var1, values_from = Rho)
PValIntronTable<-rhoAndPValIntronGenes$PValueTable %>% pivot_wider(names_from = Var1, values_from = PValue)

AllGeneRho<-rbind(RhoTable, RhoIntronTable)
AllGenePVal<-rbind(PValueTable, PValIntronTable)

write.csv(AllGeneRho, "../Output-files/CorrelationStatistics/AllGeneRho.csv", row.names = FALSE)
write.csv(AllGenePVal, "../Output-files/CorrelationStatistics/AllGenePVal.csv", row.names = FALSE)

### CDS only
setnames(ExonScores, old = names(ExonScores), new = sub("weighted_mean_NormalizedSum_|GeneMethylationRatio_", "", names(ExonScores)))

SubsetExon<-select(ExonScores, starts_with(c("H2","H3", "ATAC", "H4", "IBM", "JMJ",
                                             "Pol", "SDG", "avg", "tissue", "cv",
                                             "CG", "CHH","CHG", "first", "total","Intron")), geneLength)

epigmarks <- grep("H3|H2|H4|ATAC|^CG|^CHH|^CHG|ATAC|avg|tissue|cv|IBM|JMJ|Pol|SDG", colnames(SubsetExon), value = TRUE)

# for all genes
corAllGenes<-correlationDataTable(SubsetExon, epigmarks, c("totalintronNum", "geneLength", "totalExonLength"))

# for intron-containing genes
corIntronGenes<-correlationDataTable(SubsetExon[totalintronNum!=0], epigmarks,
                                     c("firstIntronLength", "firstintronpos",
                                       "totalIntronLength", "IntronPct"))

correlationsExons<-rbind(corAllGenes, corIntronGenes)
write.csv(correlationsExons, "../Output-files/Arabidopsis-CDSSpearmanCor.csv", row.names = FALSE)

rhoAndPValCDS<-correlationDataTableWithPValues(SubsetExon, epigmarks, c("totalintronNum", "geneLength", "totalExonLength"))
rhoTableCDS<-rhoAndPValCDS$RhoTable %>% pivot_wider(names_from = Var1, values_from = Rho)
PValueCDS<-rhoAndPValCDS$PValueTable %>% pivot_wider(names_from = Var1, values_from = PValue)

rhoAndPValCDSIntronGenes<-correlationDataTableWithPValues(SubsetExon[totalintronNum!=0], epigmarks,
                                                          c("firstIntronLength", "firstintronpos",
                                                            "totalIntronLength", "IntronPct"))
rhoCDSINtronTable<-rhoAndPValCDSIntronGenes$RhoTable%>% pivot_wider(names_from = Var1, values_from = Rho)
PValCDSINtrontable<-rhoAndPValCDSIntronGenes$PValueTable%>% pivot_wider(names_from = Var1, values_from = PValue)

CDSRho<-rbind(rhoTableCDS, rhoCDSINtronTable)
CDSPVal<-rbind(PValueCDS, PValCDSINtrontable)

write.csv(CDSRho, "../Output-files/CorrelationStatistics/CDSRho.csv", row.names = FALSE)
write.csv(CDSPVal, "../Output-files/CorrelationStatistics/CDSPVal.csv", row.names = FALSE)

### Non-coding only
setnames(IntronScores, old = names(IntronScores), new = sub("weighted_mean_NormalizedSum_|GeneMethylationRatio_", "", names(IntronScores)))

SubsetIntron<-select(IntronScores, starts_with(c("H2","H3", "ATAC", "H4", "IBM", "JMJ",
                                                 "Pol", "SDG", "avg", "tissue", "cv",
                                                 "CG", "CHH","CHG", "first", "total","Intron")), geneLength)

epigmarks <- grep("H3|H2|H4|ATAC|^CG|^CHH|^CHG|ATAC|avg|tissue|cv|IBM|JMJ|Pol|SDG", colnames(SubsetIntron), value = TRUE)

# for all genes
corAllGenes<-correlationDataTable(SubsetIntron, epigmarks, c("totalintronNum", "geneLength", "totalExonLength"))

# for intron-containing genes
corIntronGenes<-correlationDataTable(SubsetIntron[totalintronNum!=0], epigmarks,
                                     c("firstIntronLength", "firstintronpos",
                                       "totalIntronLength", "IntronPct"))

correlationsIntrons<-rbind(corAllGenes, corIntronGenes)
write.csv(correlationsIntrons, "../Output-files/Arabidopsis-NonCodingSpearmanCor.csv", row.names = FALSE)

rhoAndPValIntron<-correlationDataTableWithPValues(SubsetIntron, epigmarks, c("totalintronNum", "geneLength", "totalExonLength"))
rhoTableIntron<-rhoAndPValIntron$RhoTable %>% pivot_wider(names_from = Var1, values_from = Rho)
PValueIntron<-rhoAndPValIntron$PValueTable %>% pivot_wider(names_from = Var1, values_from = PValue)

rhoAndPValIntronIntronGenes<-correlationDataTableWithPValues(SubsetIntron[totalintronNum!=0], epigmarks,
                                                             c("firstIntronLength", "firstintronpos",
                                                               "totalIntronLength", "IntronPct"))
rhoIntronINtronTable<-rhoAndPValIntronIntronGenes$RhoTable%>% pivot_wider(names_from = Var1, values_from = Rho)
PValIntronINtrontable<-rhoAndPValIntronIntronGenes$PValueTable%>% pivot_wider(names_from = Var1, values_from = PValue)

CDSRho<-rbind(rhoTableIntron, rhoIntronINtronTable)
CDSPVal<-rbind(PValueIntron, PValIntronINtrontable)

write.csv(CDSRho, "../Output-files/CorrelationStatistics/NCSRho.csv", row.names = FALSE)
write.csv(CDSPVal, "../Output-files/CorrelationStatistics/NCSPVal.csv", row.names = FALSE)

# Random Forest ------------------------------------------------------

### All Gene
#intron-containing only
rankedTable<-rankColumns(MasterData[totalintronNum!=0], epigmarks)

marks <- grep("_ranked", colnames(rankedTable), value = TRUE)

predictors <- c("totalintronNum", "geneLength", "totalExonLength", "firstIntronLength", "firstintronpos",
                "totalIntronLength", "IntronPct")

# Example of One tree
# H3K4me1Model<-RFMarks(rankedTable, "H3K4me1_ranked", predictors)
# plot(H3K4me1Model[["H3K4me1_ranked"]])

models <- RFMarks(rankedTable, marks, predictors)

#saveRDS(models, file = "../Output-files/Arabidopsis-PCSD-RFModels1.rds") # save models as RDS object
# load using readRDS()
pctvar_values <- lapply(c(models), function(model) {
  pct_var = round(100 * model$rsq[length(model$rsq)], digits = 2)
  model$rsq[length(model$rsq)] # Last value is OOB MSE
  return(pct_var)
})

pctvar_df <- data.frame(
  Model = names(pctvar_values),
  pct_var = unlist(pctvar_values)
)

write.csv(pctvar_df, "../Output-files/Arabidopsis-PCSD-RFPctVar.csv", row.names = FALSE)

### CDS Only
rankedTableCDS<-rankColumns(ExonScores[totalintronNum!=0], epigmarks)

marks <- grep("_ranked", colnames(rankedTableCDS), value = TRUE)

predictors <- c("totalintronNum", "geneLength", "totalExonLength", "firstIntronLength", "firstintronpos",
                "totalIntronLength", "IntronPct")

modelsCDS <- RFMarks(rankedTableCDS, marks, predictors)

#saveRDS(modelsCDS, file = "../Output-files/Arabidopsis-PCSD-RFExonModels1.rds") # save models as RDS object
# load using readRDS()

pctvar_values <- lapply(c(modelsCDS), function(model) {
  pct_var = round(100 * model$rsq[length(model$rsq)], digits = 2)
  model$rsq[length(model$rsq)] # Last value is OOB MSE
  return(pct_var)
})

pctvar_df <- data.frame(
  Model = names(pctvar_values),
  pct_var = unlist(pctvar_values)
)

write.csv(pctvar_df, "../Output-files/Arabidopsis-PCSD-RFExonPctVar.csv", row.names = FALSE)

### NonCoding Only

rankedTableNonCoding<-rankColumns(IntronScores[totalintronNum!=0], epigmarks)

marks <- grep("_ranked", colnames(rankedTableNonCoding), value = TRUE)

predictors <- c("totalintronNum", "geneLength", "totalExonLength", "firstIntronLength", "firstintronpos",
                "totalIntronLength", "IntronPct")

modelsNonCoding <- RFMarks(rankedTableNonCoding, marks, predictors)

#saveRDS(modelsNonCoding, file = "../Output-files/Arabidopsis-PCSD-RFIntronModels1.rds") # save models as RDS object
# load using readRDS()

pctvar_values <- lapply(c(modelsNonCoding), function(model) {
  pct_var = round(100 * model$rsq[length(model$rsq)], digits = 2)
  model$rsq[length(model$rsq)] # Last value is OOB MSE
  return(pct_var)
})

pctvar_df <- data.frame(
  Model = names(pctvar_values),
  pct_var = unlist(pctvar_values)
)

write.csv(pctvar_df, "../Output-files/Arabidopsis-PCSD-RFIntronPctVar.csv", row.names = FALSE)


# 5X Cross Validation ----------------------------------------------------
# Process for H3K3me1, for example:
#H3K4me1CV<-processCrossvalidation(rankedTable, "H3K4me1_ranked", "1")

### All Gene 5X CV
set.seed(123)
results <- processCrossvalidation(rankedTable, marks, n_times = 5)

#saveRDS(results, file = "../Output-files/Arabidopsis-PCSD-CVModels.rds") # save models as RDS object
# load using readRDS()
# Access the results for H3K4me1, for example:
# H3K4me1_CV <- results$H3K4me1$CV
# H3K4me1_metrics <- results$H3K4me1$Metrics
# H3K4me1_features <- results$H3K4me1$Features

# Combine all Metrics into a single data frame and remove duplicates
allmetrics <- bind_rows(lapply(names(results), function(mark) {
  results[[mark]]$Metrics %>%  # Access Metrics for each mark
    mutate(mark = mark) %>%    # Add mark column to identify the source
    select(mark, RMSE, MAE, R2, COR)  # Keep relevant columns
})) %>%
  distinct()  # Remove duplicate rows

# Combine all Features into a single data frame and remove duplicates
allfeatures <- bind_rows(lapply(names(results), function(mark) {
  results[[mark]]$Features })) %>%
  distinct()  # Remove duplicate rows

write.csv(allmetrics, "../Output-files/Arabidopsis-PCSD-CVMetrics.csv", row.names = FALSE)
write.csv(allfeatures, "../Output-files/Arabidopsis-PCSD-CVFeatures.csv", row.names = FALSE)

# Combine all Metrics to obtain Predicted Values
allPredicted <- bind_rows(lapply(names(results), function(mark) {
  results[[mark]]$Metrics %>%  # Access Metrics for each mark
    mutate(mark = mark)
}))
write.csv(allPredicted, "../Output-files/Arabidopsis-PCSD-CVPredictedVals.csv", row.names = FALSE)

### CDS 5X CV
set.seed(123)
CDSresults <- processCrossvalidation(rankedTableCDS, marks, n_times = 5)

# Combine all Metrics into a single data frame and remove duplicates
allmetrics <- bind_rows(lapply(names(CDSresults), function(mark) {
  CDSresults[[mark]]$Metrics %>%  # Access Metrics for each mark
    mutate(mark = mark) %>%    # Add mark column to identify the source
    select(mark, RMSE, MAE, R2, COR)  # Keep relevant columns
})) %>%
  distinct()  # Remove duplicate rows

# Combine all Features into a single data frame and remove duplicates
allfeatures <- bind_rows(lapply(names(CDSresults), function(mark) {
  CDSresults[[mark]]$Features })) %>%
  distinct()  # Remove duplicate rows

write.csv(allmetrics, "../Output-files/Arabidopsis-PCSD-CDSCVMetrics.csv", row.names = FALSE)
write.csv(allfeatures, "../Output-files/Arabidopsis-PCSD-CDSCVFeatures.csv", row.names = FALSE)

# Combine all Metrics to obtain Predicted Values
allPredicted <- bind_rows(lapply(names(CDSresults), function(mark) {
  CDSresults[[mark]]$Metrics %>%  # Access Metrics for each mark
    mutate(mark = mark)
}))
write.csv(allPredicted, "../Output-files/Arabidopsis-PCSD-CDSCVPredictedVals.csv", row.names = FALSE)

### Non-Coding 5X CV
set.seed(123)
NCresults <- processCrossvalidation(rankedTableNonCoding, marks, n_times = 5)

# Combine all Metrics into a single data frame and remove duplicates
allmetrics <- bind_rows(lapply(names(NCresults), function(mark) {
  NCresults[[mark]]$Metrics %>%  # Access Metrics for each mark
    mutate(mark = mark) %>%    # Add mark column to identify the source
    select(mark, RMSE, MAE, R2, COR)  # Keep relevant columns
})) %>%
  distinct()  # Remove duplicate rows

# Combine all Features into a single data frame and remove duplicates
allfeatures <- bind_rows(lapply(names(NCresults), function(mark) {
  NCresults[[mark]]$Features })) %>%
  distinct()  # Remove duplicate rows

write.csv(allmetrics, "../Output-files/Arabidopsis-PCSD-NonCodingCVMetrics.csv", row.names = FALSE)
write.csv(allfeatures, "../Output-files/Arabidopsis-PCSD-NonCodingCVFeatures.csv", row.names = FALSE)

# Combine all Metrics to obtain Predicted Values
allPredicted <- bind_rows(lapply(names(NCresults), function(mark) {
  NCresults[[mark]]$Metrics %>%  # Access Metrics for each mark
    mutate(mark = mark)
}))
write.csv(allPredicted, "../Output-files/Arabidopsis-PCSD-NonCodingCVPredictedVals.csv", row.names = FALSE)


# TSS region --------------------------------------------------------------
TSSregion <- ArabidopsisFeatures$allGenes %>% select(gene,chr,start,stop,strand)
TSSregion <- TSSregion %>%
  mutate(
    min = if_else(strand == "+", start - 1500, stop - 1500),
    max = if_else(strand == "+", start + 1500, stop + 1500)
  )
# IF THIS IS - STRAND THE START SITE IS ACUALLY "STOP"

TSSregion<-TSSregion %>% select(gene, min, max, strand, chr)
names(TSSregion)[names(TSSregion) == "min"] <- "start"
names(TSSregion)[names(TSSregion) == "max"] <- "stop"
TSSregion<-isRepresentative(TSSregion, "../Input-files/TAIR10_representative_gene_models.txt")
setkey(TSSregion, chr, start, stop)
setkey(normalized, chr, start, stop)

overlap_results <- foverlaps(TSSregion, normalized, by.x = c("chr", "start", "stop"), type = "any")
overlap_results<-as.data.table(overlap_results)

overlap_results<-overlap_results %>%
  group_by(gene, strand) %>%
  mutate(n=row_number(),
         n = if_else(strand == '-',
                     rev(row_number()),
                     row_number()))

current_names <- names(overlap_results)
new_names <- gsub("NormalizedSum_", "", current_names) # Remove prefix
setnames(overlap_results, old = current_names, new = new_names)
all<-merge(ArabidopsisFeatures$IntronFeatures, overlap_results, by="gene")

write.csv(all, "../Output-files/Arabidopsis-PCSD-TSSregion.csv", row.names = FALSE)

# Gene body region --------------------------------------------------------

Genicregion <- ArabidopsisFeatures$allGenes %>% select(gene,chr,start,stop,strand)
Genicregion <- Genicregion %>%
  mutate(
    min = if_else(strand == "+", start - 3000, stop - 3000),
    max = if_else(strand == "+", start + 3000, stop + 3000)
  )
# IF THIS IS - STRAND THE START SITE IS ACUALLY "STOP"

Genicregion<-Genicregion %>% select(gene, min, max, strand, chr)
names(Genicregion)[names(Genicregion) == "min"] <- "start"
names(Genicregion)[names(Genicregion) == "max"] <- "stop"
Genicregion<-isRepresentative(Genicregion, "../Input-files/TAIR10_representative_gene_models.txt")
setkey(Genicregion, chr, start, stop)
setkey(normalized, chr, start, stop)

overlap_results <- foverlaps(Genicregion, normalized, by.x = c("chr", "start", "stop"), type = "any")
overlap_results<-as.data.table(overlap_results)

overlap_results<-overlap_results %>%
  group_by(gene, strand) %>%
  mutate(n=row_number(),
         n = if_else(strand == '-',
                     rev(row_number()),
                     row_number()))

current_names <- names(overlap_results)
new_names <- gsub("NormalizedSum_", "", current_names) # Remove prefix
setnames(overlap_results, old = current_names, new = new_names)
all<-merge(ArabidopsisFeatures$IntronFeatures, overlap_results, by="gene")

write.csv(all, "../Output-files/Arabidopsis-PCSD-GeneBodyRegion.csv", row.names = FALSE)

# Gene Duplicate Interm. Files --------------------------------------------
reps<-fread("../Input-files/TAIR10_representative_gene_models.txt")

genes<-fread("../Input-files/41586_2021_4269_MOESM4_ESM.txt")
genes$model<-gsub("ID=(.+);Par.+","\\1",genes$info)
genes$model<-gsub("Parent=(.+)","\\1",genes$model)
genes$model<-gsub("(.+),.+","\\1",genes$model)
genes$model<-gsub("-Protein","",genes$model)

genes<-genes[model %in% reps$Model,c("gene","type","chr","start","stop")]
fwrite(genes, "../Output-files/gene_features.csv")

CDS<-unique(genes)[type=="CDS"]
intron<-unique(genes)[type=="intron"]

CDS_genes<-CDS[,.(cds_length=sum(stop-start), CDS_N=.N), by=gene]
intron_genes<-intron[,.(intron_length=sum(stop-start), intron_N=.N), by=gene]

merge<-merge(CDS_genes, intron_genes, by="gene", all.x=T)
merge$intron_length[is.na(merge$intron_length)]<-0
merge$intron_N[is.na(merge$intron_N)]<-0
fwrite(merge, "../Output-files/gene_intron_CDS.csv")
gene_features<-fread("../Output-files/gene_intron_CDS.csv")

#https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.19227
dups<-fread("../Input-files/geneduplicates", skip = 1) #change this to local file
gene_features<-fread("../Output-files/gene_features.csv")
gene_intron_CDS<-fread("../Output-files/gene_intron_CDS.csv")


dups<-data.table(dup1=dups$V2, dup2=dups$V4, src=dups$V6)
x<-unlist(dups[1])
dups<-rbindlist(apply(dups, 1, function(x){
  pair<-sample(x[1:2])
  src<-x[3]
  data.table(dup1=pair[1], dup2=pair[2], src=src)
}))

dups$intron.1<-gene_intron_CDS$intron_N[match(dups$dup1, gene_intron_CDS$gene)]
dups$intronL.1<-gene_intron_CDS$intron_length[match(dups$dup1, gene_intron_CDS$gene)]

dups$intron.2<-gene_intron_CDS$intron_N[match(dups$dup2, gene_intron_CDS$gene)]
dups$intronL.2<-gene_intron_CDS$intron_length[match(dups$dup2, gene_intron_CDS$gene)]

dups$intron_diff<-dups$intron.1-dups$intron.2
dups$intron_diffL<-dups$intronL.1-dups$intronL.2

dups$lossL<-ifelse(dups$intron_diffL<0,"loss","gain")
dups$lossL[dups$intron_diffL==0]<-"equal"
dups$loss<-ifelse(dups$intron_diff<0,"loss","gain")
dups$loss[dups$intron_diff==0]<-"equal"
dups<-dups[!is.na(loss)]
fwrite(dups, "../Output-files/dups.csv")
# Gene Duplicates ---------------------------------------------------------

Dups<-fread("../Output-files/dups.csv")

current_names <- names(MasterData)
new_names <- gsub("NormalizedSum_", "", current_names) # Remove prefix
setnames(MasterData, old = current_names, new = new_names)

epigmarks <- grep("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC|avg|tissue|cv|IBM|JMJ|Pol|SDG", colnames(MasterData), value = TRUE)

all<-rbindlist(lapply(epigmarks, function(m){

  result<-intron_marks_dups(Dups, MasterData, m)
  ttest<-result$ttest
  return(data.table(mark=m, statistic=ttest$statistic, p=ttest$p.value))
}))
write.csv(all, "../Output-files/Arabidopsis-GeneDupsTTest.csv", row.names = FALSE)

# Transposed Duplicates ---------------------------------------------------
# Data from: https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.19227
Transposed<-fread("../Input-files/transposedduplicates", header = TRUE) #change this to local file
Transposed<-Transposed[,-1]

representative<-isRepresentative(ArabidopsisFeatures$IntronFeatures, "../Input-files/TAIR10_representative_gene_models.txt")
representative$gene<-gsub("\\..*","", representative$gene)

TransposedFeatures <- Transposed %>%
  inner_join(representative %>% select(gene, totalintronNum, totalIntronLength), by = c("Transposed" = "gene"))

names(TransposedFeatures)[names(TransposedFeatures) == "totalintronNum"] <- "TransposedTotalIntronNum"
names(TransposedFeatures)[names(TransposedFeatures) == "totalIntronLength"] <- "TransposedTotalIntronLength"

result2 <- Transposed %>%
  inner_join(representative %>% select(gene, totalintronNum, totalIntronLength), by = c("Parental" = "gene"))

names(result2)[names(result2) == "totalintronNum"] <- "ParentalTotalIntronNum"
names(result2)[names(result2) == "totalIntronLength"] <- "ParentalTotalIntronLength"

Merged <- merge(TransposedFeatures, result2, by = c("Transposed", "Parental", "Transposed_met", "Parental_met"))
Merged$TransposedTotalIntronLength[is.na(Merged$TransposedTotalIntronLength)] <- 0
Merged$ParentalTotalIntronLength[is.na(Merged$ParentalTotalIntronLength)] <- 0

Merged$intron_diff<-Merged$TransposedTotalIntronNum-Merged$ParentalTotalIntronNum
Merged$intron_diffL<-Merged$TransposedTotalIntronLength-Merged$ParentalTotalIntronLength

Merged$loss<-ifelse(Merged$intron_diff<0, "loss", "gain")
Merged[intron_diff==0]$loss<-"equal"

Merged$lossL<-ifelse(Merged$intron_diffL<0, "loss", "gain")
Merged[intron_diffL==0]$lossL<-"equal"
write.csv(Merged, "../Output-files/ArabidopsisTransposedFeatures.csv", row.names = F)

Tall<-rbindlist(lapply(epigmarks, function(m){

  result<-intron_marks_transposed(Merged, MasterData, m)
  ttest<-result$ttest
  return(data.table(mark=m, statistic=ttest$statistic, p=ttest$p.value))
}))
write.csv(Tall, "../Output-files/Arabidopsis-TransposedTTest.csv", row.names = FALSE)
