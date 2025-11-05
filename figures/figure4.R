# -------------------------------------------------------------------------
# File: figure4
# Pierce et al 2025
# Purpose: Create figure 4 & associated supplementary figs
# -------------------------------------------------------------------------

# Libraries ---------------------------------------------------------------
source("../functions_and_intermediates/functions.R")

# transposed dups ---------------------------------------------------------
Transposed<-fread("../Output-files/ArabidopsisTransposedFeatures.csv")

MasterData<-fread("../Output-files/Arabidopsis-PCSD-MasterData.csv")
current_names <- names(MasterData)
new_names <- gsub("NormalizedSum_", "", current_names) # Remove prefix
setnames(MasterData, old = current_names, new = new_names)

all<-fread("../Output-files/Arabidopsis-TransposedTTest.csv")

epigmarks <- grep("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC|avg|tissue|cv|IBM|JMJ|Pol|SDG", colnames(MasterData), value = TRUE)

all$mark<-factor(all$mark, levels=all$mark[rev(order(all$statistic))])

col2<-ggplot(all, aes(x=-statistic, y=mark, fill=-log10(p)))+
  geom_bar(stat="identity", colour="black")+
  theme_classic(base_size = 6)+
  theme(#axis.text.x=element_text(angle=45, hjust=1),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_gradient(low = "#edf8fb",high = "#cc4c02",
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  labs(x="T-test statistic")

plots<-lapply(epigmarks, function(m){
  intron_marks_transposed(Transposed, MasterData, m)
})

col1<-plot_grid(plots[[13]]$plot,
                plots[[3]]$plot,
                plots[[30]]$plot,
                ncol = 1, align = "v", labels = "C", label_size = 8)

col3<-plot_grid(plots[[5]]$plot,
                plots[[15]]$plot,
                plots[[4]]$plot,
                ncol = 1, align = "v", labels = "A", label_size = 8)

pdf("../Figures/fig4-rightCol.pdf", height = 3, width = 1.2)
col1
dev.off()

pdf("../Figures/fig4-leftCol.pdf", height = 3, width = 1.2)
col3
dev.off()

pdf("../Figures/fig4-ttest.pdf", height = 3.6, width = 2.2)
col2
dev.off()

col4<-plot_grid(plots[[28]]$plot,
                plots[[29]]$plot,
                plots[[27]]$plot,
                plots[[1]]$plot,ncol=1, align = "v", labels="D",
                label_size = 8)

pdf("../Figures/fig4-expnCol.pdf", height = 4, width = 1.2)
col4
dev.off()

# supp --------------------------------------------------------------------
Transposed<-fread("../Output-files/ArabidopsisTransposedFeatures.csv")

MasterData<-fread("../Output-files/Arabidopsis-PCSD-MasterData.csv")
current_names <- names(MasterData)
new_names <- gsub("NormalizedSum_", "", current_names) # Remove prefix
setnames(MasterData, old = current_names, new = new_names)

all<-fread("../Output-files/Arabidopsis-TransposedTTest.csv")

epigmarks <- grep("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC|avg|tissue|cv|IBM|JMJ|Pol|SDG", colnames(MasterData), value = TRUE)

plots<-lapply(epigmarks, function(m){
  intron_marks_transposed(Transposed, MasterData, m)
})

#exon
ExonMasterData<-fread("../Output-files/Arabidopsis-PCSD-ExonMasterData.csv")
colnames(ExonMasterData) <- gsub("weighted_mean_NormalizedSum_", "", colnames(ExonMasterData))
colnames(ExonMasterData) <- gsub("GeneMethylationRatio_", "", colnames(ExonMasterData))

plotsExon<-lapply(epigmarks, function(m){
  intron_marks_transposed(Transposed, ExonMasterData, m)
})

#intron
IntronMasterData<-fread("../Output-files/Arabidopsis-PCSD-IntronMasterData.csv")
colnames(IntronMasterData) <- gsub("weighted_mean_NormalizedSum_", "", colnames(IntronMasterData))
colnames(IntronMasterData) <- gsub("GeneMethylationRatio_", "", colnames(IntronMasterData))

plotsIntron<-lapply(epigmarks, function(m){
  intron_marks_transposed(Transposed, IntronMasterData, m)
})

#ttest
supcol1<-plot_ttestTransposed(Transposed, MasterData, epigmarks)
supcol2<-plot_ttestTransposed(Transposed, ExonMasterData, epigmarks)
supcol3<-plot_ttestTransposed(Transposed, IntronMasterData, epigmarks)

TransposedTTst<-plot_grid(supcol2, supcol3, ncol = 2, align = "h",
                          label_size = 8, labels = c("A", "B"))

pdf("../Figures/Supp4_1.pdf", height = 8, width = 6)
TransposedTTst
dev.off()

pdf("../Figures/Supp4.pdf", height = 10, width = 8)
#TransposedTTst
plot_grid(plotlist = lapply(plots, `[[`, "plot"), ncol = 4, align = "h")
plot_grid(plotlist = lapply(plotsExon, `[[`, "plot"), ncol = 4, align = "h")
plot_grid(plotlist = lapply(plotsIntron, `[[`, "plot"), ncol = 4, align = "h")
dev.off()
