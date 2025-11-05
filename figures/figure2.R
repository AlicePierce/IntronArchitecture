# -------------------------------------------------------------------------
# File: figure2
# Pierce et al 2025
# Purpose: Create figure 2 & associated supplementary figs
# -------------------------------------------------------------------------

# Libraries ---------------------------------------------------------------
source("../functions_and_intermediates/functions.R")

# Intron-Exon Plots -------------------------------------------------------

IntronExonScores<-fread("../Output-files/Arabidopsis-PCSD-IntronExonDepth.csv")
genes<-gffGenes("../Genomes/TAIR10_GFF3_genes.gff")
genes<-isRepresentative(genes, "../Input-files/TAIR10_representative_gene_models.txt")
IntronExonScores<-IntronExonScores[gene %in% genes$gene]

IntronExonScores$NumGrp<-cut2(IntronExonScores$Num, g=10)
Summary<-summarizeMarksbyGroup(IntronExonScores, grep("NormalizedSum",colnames(IntronExonScores),value=TRUE), "NumGrp")
Summary[,featureNum := paste(Feature, as.character(NumGrp), sep=" ")]
Summary$Mark<-gsub("NormalizedSum_", "", Summary$Mark)

Summary<-Summary[Mark %in% grep("H3K36me3|^H3K4", Summary$Mark,value=TRUE)]

Levels<-Summary[Mark==Mark[1]]$featureNum

Summary$featureNum<-factor(Summary$featureNum, levels = Levels)

split_data <- split(Summary, Summary$Mark)
IEPlots<-plotSplitData(split_data)
IEPlots <- lapply(seq_along(IEPlots), function(i) {
  if (i <= length(IEPlots) - 1) {
    # Remove x-axis text and title for all but the last 4 plots
    IEPlots[[i]] + theme(axis.text.x = element_blank(),
                         axis.title.x = element_blank())
  } else {
    # Keep x-axis for the last 4 plots
    IEPlots[[i]]
  }
})
col1<-plot_grid(plotlist = IEPlots, ncol = 1, align = "v", rel_heights = c(1,1,1,1.3),
                labels = "A", label_size = 8)

# TSS Metaplots -----------------------------------------------------------
TSSregion<-fread("../Output-files/Arabidopsis-PCSD-TSSregion.csv")
TSSregion<-TSSregion[totalintronNum!=0]
epigmarks<-grep("H3K4|H3K36me3", colnames(TSSregion), value = TRUE)

firstintronposPlots<-featuresMetaplots(TSSregion, "firstintronpos", marks = epigmarks)
pdf("../Figures/Fig2_LegendIntronPos.pdf", height = 5.4, width = 4.5)
firstintronposPlots[[1]]
dev.off()

firstintronposPlots <- lapply(seq_along(firstintronposPlots), function(i) {
  if (i <= length(firstintronposPlots) - 1) {
    # Remove x-axis text and title for all but the last 4 plots
    firstintronposPlots[[i]] + theme(legend.position = "none",
                                     axis.text.x = element_blank())
    #axis.title.y = element_blank())
  } else {
    # Keep x-axis for the last 4 plots
    firstintronposPlots[[i]] + theme(#axis.title.y = element_blank(),
      legend.position = "none")
  }
})
col2<-plot_grid(plotlist = firstintronposPlots, ncol = 1, align = "v")


# GB Metaplots ------------------------------------------------------------
GBregion<-fread("../Output-files/Arabidopsis-PCSD-GeneBodyRegion.csv")
GBregion<-GBregion[totalintronNum!=0]
epigmarks<-grep("H3K4|H3K36me3", colnames(GBregion), value = TRUE)

intronNumPlots<-featuresGeneMetaplots(GBregion, "totalintronNum", marks = epigmarks)

pdf("../Figures/Fig2_LegendIntronNum.pdf", height = 5.4, width = 4.5)
intronNumPlots[[1]]
dev.off()

intronNumPlots <- lapply(seq_along(intronNumPlots), function(i) {
  if (i <= length(intronNumPlots) - 1) {
    # Remove x-axis text and title for all but the last 4 plots
    intronNumPlots[[i]] + theme(legend.position = "none",
                                axis.text.x = element_blank(),
                                axis.title.y = element_blank())
  } else {
    # Keep x-axis for the last 4 plots
    intronNumPlots[[i]] + theme(axis.title.y=element_blank(),
                                legend.position = "none")
  }
})
col3<-plot_grid(plotlist = intronNumPlots, ncol = 1, align = "v")

pdf("../Figures/Fig2_IntronExon.pdf", height = 6, width = 2.5)
col1
dev.off()

metaplots<-plot_grid(col2, col3, ncol=2, labels = c("B", "C"), label_size = 8, rel_widths = c(1,2))

pdf("../Figures/Fig2_Metaplots_noLegend.pdf", height = 5.4, width = 4.5)
metaplots
dev.off()

fig2<-plot_grid(col1, col2, col3, ncol=3, labels = c("A", "B", "C"),
                label_size = 8, rel_widths = c(2,2,3), align = "v", axis = "bt")

pdf("../Figures/Fig2Mockup2.pdf", height = 6, width = 7)
fig2
dev.off()

# SuppFigs ----------------------------------------------------------------

#intron exon plots
IntronExonScores<-fread("../Output-files/Arabidopsis-PCSD-IntronExonDepth.csv")
genes<-gffGenes("../Genomes/TAIR10_GFF3_genes.gff")
genes<-isRepresentative(genes, "../Input-files/TAIR10_representative_gene_models.txt")
IntronExonScores<-IntronExonScores[gene %in% genes$gene]

IntronExonScores$NumGrp<-cut2(IntronExonScores$Num, g=10)
Summary<-summarizeMarksbyGroup(IntronExonScores, grep("NormalizedSum",colnames(IntronExonScores),value=TRUE), "NumGrp")
Summary[,featureNum := paste(Feature, as.character(NumGrp), sep=" ")]

Levels<-Summary[Mark==Mark[1]]$featureNum

Summary$featureNum<-factor(Summary$featureNum, levels = Levels)
Summary$Mark<-gsub("NormalizedSum_", "", Summary$Mark)

split_data <- split(Summary, Summary$Mark)

AllIEPlots<-plotSplitData(split_data)

col1<-plot_grid(plotlist = AllIEPlots, ncol = 4)

pdf("../Figures/Supp2_IntronExonScatter.pdf", height = 10, width = 8)
col1
dev.off()
