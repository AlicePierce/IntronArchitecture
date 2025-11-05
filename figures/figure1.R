# -------------------------------------------------------------------------
# File: figure1
# Pierce et al 2025
# Purpose: Create figure 1 & associated supplementary figs
# -------------------------------------------------------------------------

# Libraries ---------------------------------------------------------------
source("../functions_and_intermediates/functions.R")

# %VarExp -----------------------------------------------------------------
PctVar<-fread("../Output-files/Arabidopsis-PCSD-RFPctVar.csv")
PctVar$Model<-gsub("_ranked", "", PctVar$Model)

PctVar <- PctVar %>%
  arrange(pct_var)

marks <- grep("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC", PctVar$Model, value = TRUE)
expn<-grep("avg|tissue|cv", PctVar$Model, value=TRUE)
prot<-grep("IBM|SDG|Pol|JMJ", PctVar$Model, value=TRUE)

PctVar$Model<-factor(PctVar$Model, levels=c(expn,prot, marks))

# Percent Variance Plot
PctVarMarks<-ggplot(PctVar, aes(x =-pct_var, y =  Model)) +
  geom_bar(stat = "identity", fill = "steelblue", color="black", width = 0.75) +
  theme_classic(base_size = 6) +
  labs(x = "% Var Explained") +
  scale_x_continuous(labels = abs, limits = c(-48,4.5))+
  scale_y_discrete(position = "right") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(hjust = 0.5),
        axis.line.y = element_blank(),
        strip.text = element_blank())

# Heatmap -----------------------------------------------------------------
Correlations<-fread("../Output-files/CorrelationStatistics/AllGeneRho.csv")
Correlations<-Correlations %>% pivot_longer(cols=!Var2,names_to = "mark", values_to = "value")
colnames(Correlations)<-c("rn", "mark", "value")

RFFeatures<-fread("../Output-files/Arabidopsis-PCSD-CVFeatures.csv")
featureImp<-RFFeatures %>% group_by(rn,mark) %>%
  summarise(`%IncMSE`=mean(`%IncMSE`),
            IncNodePurity=mean(IncNodePurity),
            n=n())
featureImp$mark<-gsub("_ranked", "", featureImp$mark)

DT<-merge(Correlations, featureImp, by=c("mark", "rn"))
DT$mark<-factor(DT$mark, levels=c(expn, prot, marks))

# Rename features
DT$rn[DT$rn == "firstIntronLength"] <- "First Intron Length"
DT$rn[DT$rn == "firstintronpos"] <- "First Intron Position"
DT$rn[DT$rn == "geneLength"] <- "Gene Length"
DT$rn[DT$rn == "IntronPct"] <- "Intron Percent"
DT$rn[DT$rn == "totalExonLength"] <- "All Exon Length"
DT$rn[DT$rn == "totalIntronLength"] <- "All Intron Length"
DT$rn[DT$rn == "totalintronNum"] <- "Number of Introns"

# Heatmap plot
CorPlot<-ggplot(DT, aes(x=rn, y=mark, fill=value))+
  geom_point(aes(size = `%IncMSE`), col="black",shape = 21, stroke = 0.5) +
  coord_fixed()+
  theme_classic(base_size = 6)+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.key.size = unit(0.25, "cm"),
        axis.text.y = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA))+ # Transparent plot background)+
  scale_fill_gradient2(high = "#0570b0", low = "#cc4c02", mid = "#fff7fb", name = "Spearman \ncoefficient",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))

# Save Percent Var and Heatmap as one figure
pdf("../Figures/CorrelationAndRF.pdf", height = 6.75, width = 4.75)
plot_grid(PctVarMarks, CorPlot, labels = c("A", "B"), label_size = 8,
          rel_heights = 1, align = "h")
dev.off()

# Scatterplots ------------------------------------------------------------
Marks<-fread("../Output-files/Arabidopsis-PCSD-MasterData.csv")
colnames(Marks) <- gsub("^NormalizedSum_", "", colnames(Marks))

epigmarks<-grep("H3K4me1|H3K4me3|H3K36me3|H3K27me3", colnames(Marks), value = TRUE)
order<-c("H3K4me3", "H3K4me1", "H3K36me3", "H3K27me3")
epigmarks<-factor(epigmarks, levels = order)

IntronContGenes<-Marks[totalintronNum!=0]

# Intron Number
Marks$intronGroup <- cut2(Marks$totalintronNum, g=10)
all<-introngroup_summaries(Marks, epigmarks)
all$Mark<-factor(all$Mark, levels=epigmarks)
plots_list1<-lapply(order, function(mark) plot_introngroup2(all, mark, "Number of Introns"))

plots_list1 <- lapply(seq_along(plots_list1), function(i) {
  if (i <= length(plots_list1) - 1) {
    plots_list1[[i]] + theme(axis.text.x = element_blank(),
                             axis.title.x = element_blank())
  } else {
    plots_list1[[i]]+
      theme(axis.text.x=element_text(angle=45, hjust=1))
  }
})

column3 <- plot_grid(plotlist = plots_list1, ncol = 1, align = "v", rel_heights = c(1,1,1,1.2))

# First Intron Pos
IntronContGenes$intronGroup <- cut2(IntronContGenes$firstintronpos, g=10)
all<-introngroup_summaries(IntronContGenes, epigmarks)
all$Mark<-factor(all$Mark, levels=epigmarks)
plots_list2<-lapply(order, function(mark) plot_introngroup2(all, mark, "First Intron Position"))

plots_list2 <- lapply(seq_along(plots_list2), function(i) {
  if (i <= length(plots_list2) - 1) {
    plots_list2[[i]] + theme(axis.text.x = element_blank(),
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank())
  } else {
    plots_list2[[i]]+
      theme(axis.text.x=element_text(angle=45, hjust=1),
            axis.title.y = element_blank())
  }
})

column4 <- plot_grid(plotlist=plots_list2, ncol = 1, align = "v",
                     rel_heights = c(1,1,1,1.2))

scatterplots<-plot_grid(column3, column4, align = "v", axis="l",
                        labels = c("C", "D"), label_size = 8)

pdf("../Figures/Fig1-ScatterPlots.pdf", height = 5.7, width = 3)
scatterplots
dev.off()

# Supp --------------------------------------------------------------------

# Pct Var
CDSPctVar<-fread("../Output-files/Arabidopsis-PCSD-RFExonPctVar.csv")
CDSPctVarPlot<-pctVarPlot(CDSPctVar)

CDSPctVar$Model<-gsub("_ranked", "", CDSPctVar$Model)

CDSmarks <- CDSPctVar$Model[grepl("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC", CDSPctVar$Model)][order(CDSPctVar$pct_var[grepl("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC", CDSPctVar$Model)])]
CDSexpn <- CDSPctVar$Model[grepl("avg|tissue|cv", CDSPctVar$Model)][order(CDSPctVar$pct_var[grepl("avg|tissue|cv", CDSPctVar$Model)])]
CDSproteins <- CDSPctVar$Model[grepl("IBM|SDG|Pol|JMJ", CDSPctVar$Model)][order(CDSPctVar$pct_var[grepl("IBM|SDG|Pol|JMJ", CDSPctVar$Model)])]

NCPctVar<-fread("../Output-files/Arabidopsis-PCSD-RFIntronPctVar.csv")
NCPctVarPlot<-pctVarPlot(NCPctVar)

NCPctVar$Model<-gsub("_ranked", "", NCPctVar$Model)

NCmarks <- NCPctVar$Model[grepl("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC", NCPctVar$Model)][order(NCPctVar$pct_var[grepl("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC", NCPctVar$Model)])]
NCexpn <- NCPctVar$Model[grepl("avg|tissue|cv", NCPctVar$Model)][order(NCPctVar$pct_var[grepl("avg|tissue|cv", NCPctVar$Model)])]
NCproteins <- NCPctVar$Model[grepl("IBM|SDG|Pol|JMJ", NCPctVar$Model)][order(NCPctVar$pct_var[grepl("IBM|SDG|Pol|JMJ", NCPctVar$Model)])]

# Heatmaps
## CDS
CorrelationsExons<-fread("../Output-files/CorrelationStatistics/CDSRho.csv")
CorrelationsExons<-CorrelationsExons %>% pivot_longer(cols=!Var2,names_to = "mark", values_to = "value")
colnames(CorrelationsExons)<-c("rn", "mark", "value")

#CorrelationsExons$rn<-gsub("_ranked", "", CorrelationsExons$rn)

CDSRFFeatures<-fread("../Output-files/Arabidopsis-PCSD-CDSCVFeatures.csv")
CDSfeatureImp<-CDSRFFeatures %>% group_by(rn,mark) %>%
  summarise(`%IncMSE`=mean(`%IncMSE`),
            IncNodePurity=mean(IncNodePurity),
            n=n())
CDSfeatureImp$mark<-gsub("_ranked", "", CDSfeatureImp$mark)

ExonDT<-merge(CorrelationsExons, CDSfeatureImp, by=c("mark", "rn"))
ExonDT$mark<-factor(ExonDT$mark, levels=c(CDSexpn, CDSproteins, CDSmarks))

ExonDT$rn[ExonDT$rn == "firstIntronLength"] <- "First Intron Length"
ExonDT$rn[ExonDT$rn == "firstintronpos"] <- "First Intron Position"
ExonDT$rn[ExonDT$rn == "geneLength"] <- "Gene Length"
ExonDT$rn[ExonDT$rn == "IntronPct"] <- "Intron Percent"
ExonDT$rn[ExonDT$rn == "totalExonLength"] <- "All Exon Length"
ExonDT$rn[ExonDT$rn == "totalIntronLength"] <- "All Intron Length"
ExonDT$rn[ExonDT$rn == "totalintronNum"] <- "Number of Introns"

CDSMap<-ggplot(ExonDT, aes(x=rn, y=mark, fill=value))+
  geom_point(aes(size = `%IncMSE`), col="black",shape = 21, stroke = 0.5) +
  coord_fixed()+
  theme_classic(base_size = 6)+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.key.size = unit(0.25, "cm"),
        axis.text.y = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_gradient2(high = "#0570b0", low = "#cc4c02", mid = "#fff7fb",limits = c(-0.6, 0.8) , name = "Spearman \ncoefficient",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  scale_size_continuous(range = c(2, 6), limits = c(10, 90))

## Non-Coding
CorrelationsIntrons<-fread("../Output-files/CorrelationStatistics/NCSRho.csv")
CorrelationsIntrons<-CorrelationsIntrons %>% pivot_longer(cols=!Var2,names_to = "mark", values_to = "value")
colnames(CorrelationsIntrons)<-c("rn", "mark", "value")

NonCodRFFeatures<-fread("../Output-files/Arabidopsis-PCSD-NonCodingCVFeatures.csv")
NCfeatureImp<-NonCodRFFeatures %>% group_by(rn,mark) %>%
  summarise(`%IncMSE`=mean(`%IncMSE`),
            IncNodePurity=mean(IncNodePurity),
            n=n())
NCfeatureImp$mark<-gsub("_ranked", "", NCfeatureImp$mark)

IntronDT<-merge(CorrelationsIntrons, NCfeatureImp, by=c("mark", "rn"))

IntronDT$mark<-factor(IntronDT$mark, levels=c(NCexpn, NCproteins, NCmarks))

IntronDT$rn[IntronDT$rn == "firstIntronLength"] <- "First Intron Length"
IntronDT$rn[IntronDT$rn == "firstintronpos"] <- "First Intron Position"
IntronDT$rn[IntronDT$rn == "geneLength"] <- "Gene Length"
IntronDT$rn[IntronDT$rn == "IntronPct"] <- "Intron Percent"
IntronDT$rn[IntronDT$rn == "totalExonLength"] <- "All Exon Length"
IntronDT$rn[IntronDT$rn == "totalIntronLength"] <- "All Intron Length"
IntronDT$rn[IntronDT$rn == "totalintronNum"] <- "Number of Introns"

NCMap<-ggplot(IntronDT, aes(x=rn, y=mark, fill=value))+
  geom_point(aes(size = `%IncMSE`), col="black",shape = 21, stroke = 0.5) +
  coord_fixed()+
  theme_classic(base_size = 6)+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.key.size = unit(0.25, "cm"),
        axis.text.y = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_gradient2(high = "#0570b0", low = "#cc4c02", mid = "#fff7fb",limits = c(-0.6, 0.8) , name = "Spearman \ncoefficient",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  scale_size_continuous(range = c(2, 6), limits = c(10, 90))

row1<-plot_grid(CDSPctVarPlot, CDSMap, nrow = 1, align="h", rel_heights = c(1,1),axis = "l")
row2<-plot_grid(NCPctVarPlot, NCMap, nrow = 1, align = "h", rel_heights = c(1,1),axis = "l")

pdf("../Figures/Supp1_Heatmaps.pdf", height = 6.5, width = 4)
row1
row2
dev.off()

# Actual vs Predicted scatterplots # Keeping for QC
# Predicted<-fread("../Output-files/Arabidopsis-PCSD-CVPredictedVals.csv")
# ggplot(Predicted, aes(x=))

# ggplot(H3K14ac_CV[[1]][1][[1]], aes(x=H3K14ac_ranked, y=predicted))+
#   geom_point()+
#   theme_classic()

# RF CV Metrics
allMetrics<-fread("../Output-files/Arabidopsis-PCSD-CVMetrics.csv")
allMetrics$mark<-gsub("_ranked", "", allMetrics$mark)
allMetrics$mark <- gsub("tissue_breadth", "Tissue Breadth", allMetrics$mark)
allMetrics$mark<- gsub("avgExpression", "Avg. Expression", allMetrics$mark)
allMetrics$mark <- gsub("cv", "CV", allMetrics$mark)

plotsAll<-list(
  metricsViolinPlot(allMetrics, "RMSE"),
  metricsViolinPlot(allMetrics, "MAE"),
  metricsViolinPlot(allMetrics, "R2"),
  metricsViolinPlot(allMetrics, "COR"))

exonMetrics<-fread("../Output-files/Arabidopsis-PCSD-CDSCVMetrics.csv")
exonMetrics$mark<-gsub("_ranked", "", exonMetrics$mark)
exonMetrics$mark<- gsub("tissue_breadth", "Tissue Breadth", exonMetrics$mark)
exonMetrics$mark<- gsub("avgExpression", "Avg. Expression", exonMetrics$mark)
exonMetrics$mark <- gsub("cv", "CV", exonMetrics$mark)

plotsExon<-list(
  metricsViolinPlot(exonMetrics, "RMSE"),
  metricsViolinPlot(exonMetrics, "MAE"),
  metricsViolinPlot(exonMetrics, "R2"),
  metricsViolinPlot(exonMetrics, "COR"))

intronMetrics<-fread("../Output-files/Arabidopsis-PCSD-NonCodingCVMetrics.csv")
intronMetrics$mark<-gsub("_ranked", "", intronMetrics$mark)
intronMetrics$mark<- gsub("tissue_breadth", "Tissue Breadth", intronMetrics$mark)
intronMetrics$mark<- gsub("avgExpression", "Avg. Expression", intronMetrics$mark)
intronMetrics$mark <- gsub("cv", "CV", intronMetrics$mark)

plotsIntron<-list(
  metricsViolinPlot(intronMetrics, "RMSE"),
  metricsViolinPlot(intronMetrics, "MAE"),
  metricsViolinPlot(intronMetrics, "R2"),
  metricsViolinPlot(intronMetrics, "COR"))

metricsAllgrid<-plot_grid(plotlist = plotsAll, ncol = 2, align = "v")
metricsExongrid<-plot_grid(plotlist = plotsExon, ncol = 2, align = "v")
metricsIntrongrid<-plot_grid(plotlist = plotsIntron, ncol = 2, align = "v")

pdf("../Figures/RFCV_supps.pdf", height = 3,width = 7)
metricsAllgrid
metricsExongrid
metricsIntrongrid
dev.off()

