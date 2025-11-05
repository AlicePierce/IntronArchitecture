# -------------------------------------------------------------------------
# File: scatter-shiny-intermediates
# Pierce et al 2025
# Purpose: Create intermediate files for scatterplots_shiny
# -------------------------------------------------------------------------

# Libraries ---------------------------------------------------------------
source("functions.R")

# Scatter Data Frame ------------------------------------------------------

### All Gene
Marks<-fread("../Output-files/Arabidopsis-PCSD-MasterData.csv")
colnames(Marks) <- gsub("^NormalizedSum_", "", colnames(Marks))
colnames(Marks) <- gsub("tissue_breadth", "Tissue Breadth", colnames(Marks))
colnames(Marks) <- gsub("avgExpression", "Avg. Expression", colnames(Marks))
colnames(Marks) <- gsub("cv", "CV", colnames(Marks))

epigmarks<-grep("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC|Avg|Tissue|CV|IBM|SDG|Pol|JMJ", colnames(Marks), value = TRUE)

IntronContGenes<-Marks[totalintronNum!=0]

gp1<-IntronGroupDF(Marks, "totalintronNum", epigmarks)
gp2<-IntronGroupDF(IntronContGenes, "firstintronpos", epigmarks)
gp3<-IntronGroupDF(IntronContGenes, "IntronPct", epigmarks)
gp4<-IntronGroupDF(IntronContGenes, "firstIntronLength", epigmarks)
gp5<-IntronGroupDF(Marks, "geneLength", epigmarks)
gp6<-IntronGroupDF(Marks, "totalExonLength", epigmarks)
gp7<-IntronGroupDF(IntronContGenes, "totalIntronLength", epigmarks)

all<-rbind(gp1, gp2, gp3, gp4, gp5, gp6, gp7)
all$region<-"AllGene"

### Exon
MarksCDS<-fread("../../Output-files/Arabidopsis-PCSD-ExonMasterData.csv")
colnames(MarksCDS) <- gsub("weighted_mean_NormalizedSum_", "", colnames(MarksCDS))
colnames(MarksCDS) <- gsub("GeneMethylationRatio_", "", colnames(MarksCDS))
colnames(MarksCDS) <- gsub("tissue_breadth", "Tissue Breadth", colnames(MarksCDS))
colnames(MarksCDS) <- gsub("avgExpression", "Avg. Expression", colnames(MarksCDS))
colnames(MarksCDS) <- gsub("cv", "CV", colnames(MarksCDS))

epigmarks<-grep("H3|H2|H4|ATAC|^CG|^CHH|^CHG|ATAC|Avg|Tissue|CV|IBM|SDG|Pol|JMJ", colnames(MarksCDS), value = TRUE)

IntronContGenesCDS<-MarksCDS[totalintronNum!=0]

gp8<-IntronGroupDF(MarksCDS, "totalintronNum", epigmarks)
gp9<-IntronGroupDF(IntronContGenesCDS, "firstintronpos", epigmarks)
gp10<-IntronGroupDF(IntronContGenesCDS, "IntronPct", epigmarks)
gp11<-IntronGroupDF(IntronContGenesCDS, "firstIntronLength", epigmarks)
gp12<-IntronGroupDF(MarksCDS, "geneLength", epigmarks)
gp13<-IntronGroupDF(MarksCDS, "totalExonLength", epigmarks)
gp14<-IntronGroupDF(IntronContGenesCDS, "totalIntronLength", epigmarks)

allExon<-rbind(gp8, gp9, gp10, gp11, gp12, gp13, gp14)
allExon$region<-"Exon"

### Intron
MarksNCS<-fread("../../Output-files/Arabidopsis-PCSD-IntronMasterData.csv") #only intron-containing genes
colnames(MarksNCS) <- gsub("weighted_mean_NormalizedSum_", "", colnames(MarksNCS))
colnames(MarksNCS) <- gsub("GeneMethylationRatio_", "", colnames(MarksNCS))
colnames(MarksNCS) <- gsub("tissue_breadth", "Tissue Breadth", colnames(MarksNCS))
colnames(MarksNCS) <- gsub("avgExpression", "Avg. Expression", colnames(MarksNCS))
colnames(MarksNCS) <- gsub("cv", "CV", colnames(MarksNCS))

epigmarks<-grep("H3|H2|H4|ATAC|^CG|^CHH|^CHG|ATAC|Avg|Tissue|CV|IBM|SDG|Pol|JMJ", colnames(MarksNCS), value = TRUE)

gp15<-IntronGroupDF(MarksNCS, "totalintronNum", epigmarks)
gp16<-IntronGroupDF(MarksNCS, "firstintronpos", epigmarks)
gp17<-IntronGroupDF(MarksNCS, "IntronPct", epigmarks)
gp18<-IntronGroupDF(MarksNCS, "firstIntronLength", epigmarks)
gp19<-IntronGroupDF(MarksNCS, "geneLength", epigmarks)
gp20<-IntronGroupDF(MarksNCS, "totalExonLength", epigmarks)
gp21<-IntronGroupDF(MarksNCS, "totalIntronLength", epigmarks)

allIntron<-rbind(gp15, gp16, gp17, gp18, gp19, gp20, gp21)
allIntron$region<-"Intron"

ScatterData<-rbind(all, allExon, allIntron)
saveRDS(ScatterData, "ScatterData.rds")
