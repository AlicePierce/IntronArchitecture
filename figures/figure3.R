# -------------------------------------------------------------------------
# File: figure3
# Pierce et al 2025
# Purpose: Create figure 3 & associated supplementary figs
# -------------------------------------------------------------------------

# Libraries ---------------------------------------------------------------
source("../functions_and_intermediates/functions.R")

# Gene duplicates ---------------------------------------------------------
Dups<-fread("../Input-files/dups.csv")
MasterData<-fread("../Output-files/Arabidopsis-PCSD-MasterData.csv")
current_names <- names(MasterData)
new_names <- gsub("NormalizedSum_", "", current_names) # Remove prefix
setnames(MasterData, old = current_names, new = new_names)
all<-fread("../Output-files/Arabidopsis-GeneDupsTTest.csv")

epigmarks <- grep("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC|avg|tissue|cv|IBM|JMJ|Pol|SDG", colnames(MasterData), value = TRUE)

all$mark<-factor(all$mark, levels=all$mark[order(all$statistic)])

ttest<-ggplot(all, aes(x=mark, y=-statistic, fill=-log10(p)))+
  geom_bar(stat="identity", colour="black")+
  theme_classic(base_size = 6)+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_gradient(low = "#edf8fb",
                      high = "#0570b0",
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  labs(y="T-test statistic")

plots<-lapply(epigmarks, function(m){
  intron_marks_dups(Dups, MasterData, m)
})

row1<-plot_grid(plots[[13]]$plot,
                plots[[3]]$plot,
                plots[[30]]$plot,
                plots[[14]]$plot,nrow = 1, align = "v", labels = "A",
                label_size = 8)

row2<-plot_grid(plots[[11]]$plot,
                plots[[15]]$plot,
                plots[[5]]$plot,
                plots[[4]]$plot,nrow = 1, align = "v", labels = "C",
                label_size = 8)

pdf("../Figures/fig3-topRow.pdf", height = 1, width = 4.5)
row1
dev.off()

pdf("../Figures/fig3-bottomRow.pdf", height = 1, width = 4.5)
row2
dev.off()

pdf("../Figures/fig3-ttest.pdf", height = 1.8, width = 3.8)
ttest
dev.off()

col2<-plot_grid(plots[[28]]$plot,
                plots[[29]]$plot,
                plots[[27]]$plot,
                plots[[1]]$plot,ncol=1, align = "v", labels = "D",
                label_size = 8)

pdf("../Figures/fig3-rightCol.pdf", height = 4, width = 1.2)
col2
dev.off()

# supp --------------------------------------------------------------------
Dups<-fread("../Input-files/dups.csv")
MasterData<-fread("../Output-files/Arabidopsis-PCSD-MasterData.csv")
current_names <- names(MasterData)
new_names <- gsub("NormalizedSum_", "", current_names) # Remove prefix
setnames(MasterData, old = current_names, new = new_names)
all<-fread("../Output-files/Arabidopsis-GeneDupsTTest.csv")

epigmarks <- grep("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC|avg|tissue|cv|IBM|JMJ|Pol|SDG", colnames(MasterData), value = TRUE)

plots<-lapply(epigmarks, function(m){
  intron_marks_dups(Dups, MasterData, m)
})

#ttest
wgd<-plot_ttest_by_source("wgd", Dups, MasterData, epigmarks)
proximal<-plot_ttest_by_source("proximal", Dups, MasterData, epigmarks)
dispersed<-plot_ttest_by_source("dispersed", Dups, MasterData, epigmarks)
tandem<-plot_ttest_by_source("tandem", Dups, MasterData, epigmarks)
transposed<-plot_ttest_by_source("transposed", Dups, MasterData, epigmarks)

dupsrc<-plot_grid(wgd, proximal, dispersed, tandem, transposed,
                  nrow = 5, align = "h", label_size = 8,
                  labels = c("A", "B", "C", "D", "E"))

#ttest exon
ExonMasterData<-fread("../Output-files/Arabidopsis-PCSD-ExonMasterData.csv")
colnames(ExonMasterData) <- gsub("weighted_mean_NormalizedSum_", "", colnames(ExonMasterData))
colnames(ExonMasterData) <- gsub("GeneMethylationRatio_", "", colnames(ExonMasterData))

wgdExon<-plot_ttest_by_source("wgd", Dups, ExonMasterData, epigmarks)
proximalExon<-plot_ttest_by_source("proximal", Dups, ExonMasterData, epigmarks)
dispersedExon<-plot_ttest_by_source("dispersed", Dups, ExonMasterData, epigmarks)
tandemExon<-plot_ttest_by_source("tandem", Dups, ExonMasterData, epigmarks)
transposedExon<-plot_ttest_by_source("transposed", Dups, ExonMasterData, epigmarks)

dupsrcExon<-plot_grid(wgdExon, proximalExon, dispersedExon, tandemExon, transposedExon,
                      nrow = 5, align = "h", label_size = 8,
                      labels = c("A", "B", "C", "D", "E"))

plotsExon<-lapply(epigmarks, function(m){
  intron_marks_dups(Dups, ExonMasterData, m)
})
#ttest intron
IntronMasterData<-fread("../Output-files/Arabidopsis-PCSD-IntronMasterData.csv")
colnames(IntronMasterData) <- gsub("weighted_mean_NormalizedSum_", "", colnames(IntronMasterData))
colnames(IntronMasterData) <- gsub("GeneMethylationRatio_", "", colnames(IntronMasterData))

wgdIntron<-plot_ttest_by_source("wgd", Dups, IntronMasterData, epigmarks)
proximalIntron<-plot_ttest_by_source("proximal", Dups, IntronMasterData, epigmarks)
dispersedIntron<-plot_ttest_by_source("dispersed", Dups, IntronMasterData, epigmarks)
tandemIntron<-plot_ttest_by_source("tandem", Dups, IntronMasterData, epigmarks)
transposedIntron<-plot_ttest_by_source("transposed", Dups, IntronMasterData, epigmarks)

dupsrcIntron<-plot_grid(wgdIntron, proximalIntron, dispersedIntron, tandemIntron, transposedIntron,
                        nrow = 5, align = "h", label_size = 8,
                        labels = c("A", "B", "C", "D", "E"))

plotsIntron<-lapply(epigmarks, function(m){
  intron_marks_dups(Dups, IntronMasterData, m)
})

Sources<-plot_grid(dupsrc, dupsrcExon, dupsrcIntron,
                   ncol = 3, align = "h", label_size = 8,
                   labels = c("All Gene", "CDS", "NC"))

#ttest
supcol1<-plot_ttest(Dups, MasterData, epigmarks)
supcol2<-plot_ttest(Dups, ExonMasterData, epigmarks)
supcol3<-plot_ttest(Dups, IntronMasterData, epigmarks)

DupsTTst<-plot_grid(supcol1, supcol2, supcol3, nrow = 3, align = "h",
                    label_size = 8, labels = c("A", "B", "C"))

pdf("../Figures/Supp3.pdf", height = 10, width = 8)
DupsTTst
dupsrc
dupsrcExon
dupsrcIntron
plot_grid(plotlist = lapply(plots, `[[`, "plot"), ncol = 4, align = "h")
plot_grid(plotlist = lapply(plotsExon, `[[`, "plot"), ncol = 4, align = "h")
plot_grid(plotlist = lapply(plotsIntron, `[[`, "plot"), ncol = 4, align = "h")
dev.off()
