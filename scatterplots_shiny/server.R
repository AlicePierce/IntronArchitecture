# -------------------------------------------------------------------------
# File: server.R
# Pierce et al 2025
# Purpose: Create Scatterplot Shiny App
# -------------------------------------------------------------------------

library(shiny)
library(ggplot2)
library(bslib)

data<-readRDS("ScatterData.rds")

server <- function(input, output) {
  output$dynamicHeader <- renderText({
    paste(input$selectGeneFeature, "~", input$selectChromatinFeature)
  })

  output$scatterplot <- renderPlot({
    region <- switch(input$selectRegion,
                     "All Gene" = "AllGene",
                     "All Exon" = "Exon",
                     "All Intron" = "Intron"
    )

    feature<-switch (input$selectGeneFeature,
                     "Intron Number" = "totalintronNum",
                     "First Intron Position" = "firstintronpos",
                     "Intron Percent" = "IntronPct",
                     "First Intron Length" = "firstIntronLength",
                     "Gene Length" = "geneLength",
                     "Total Exon Length" = "totalExonLength",
                     "Total Intron Length" = "totalIntronLength"
    )

    mark<-switch(input$selectChromatinFeature,
                 "ATAC"="ATAC",
                 "Average Expression"="Avg. Expression",
                 "CG"="CG",
                 "CHG"="CHG",
                 "CHH"="CHH",
                 "CV"="CV",
                 "H2A.W"="H2A.W",
                 "H2A.X"="H2A.X",
                 "H2A.Z"="H2A.Z",
                 "H3.1"="H3.1",
                 "H3.3"="H3.3",
                 "H3K14ac"="H3K14ac",
                 "H3K23ac"="H3K23ac",
                 "H3K27ac"="H3K27ac",
                 "H3K27me1"="H3K27me1",
                 "H3K27me3"="H3K27me3",
                 "H3K36ac"="H3K36ac",
                 "H3K36me3"="H3K36me3",
                 "H3K4me1"="H3K4me1",
                 "H3K4me2"="H3K4me2",
                 "H3K4me3"="H3K4me3",
                 "H3K56ac"="H3K56ac",
                 "H3K9ac"="H3K9ac",
                 "H3K9me1"="H3K9me1",
                 "H3K9me2"="H3K9me2",
                 "H4K16ac"="H4K16ac",
                 "IBM1"="IBM1",
                 "JMJ14"="JMJ14",
                 "Pol2"="Pol2",
                 "PolIV"="PolIV",
                 "SDG8"="SDG8",
                 "Tissue Breadth"="Tissue Breadth")

    ggplot(data[data$region == region & data$Mark == mark & data$Var==feature, ], aes(x = intronGroup, y = Mean)) +
      geom_point(shape = 21, size= 4,fill="#cc4c02") +
      geom_linerange(aes(ymin = Mean - SE, ymax = Mean + SE), color = "black", alpha = 0.6) +
      theme_classic(base_size = 24) +
      labs(x = input$selectGeneFeature, y = input$selectChromatinFeature)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
}
