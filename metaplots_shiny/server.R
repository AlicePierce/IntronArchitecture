# File: server.R
# Pierce et al 2025
# Purpose: Create Metaplot Shiny App
# -------------------------------------------------------------------------

library(shiny)
library(ggplot2)
library(bslib)

data<-readRDS("metaplotsData.rds")

server <- function(input, output) {

  output$dynamicHeader <- renderText({
    paste(input$selectChromatinFeature, "binned by",input$selectGeneFeature)
  })

  output$metaplot <- renderPlot({

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
                 "SDG8"="SDG8")


    ggplot(data[data$Mark == mark & data$geneFeature==feature, ],
           aes(x = n, y = Mean, col = intronGp, group = intronGp)) +
      geom_line(linewidth = 0.8) +
      theme_classic(base_size = 16) +
      scale_color_manual(values = c("#d0d1e9", "#a6bddb", "#74a9cf", "#3690c0",
                                    "#0570b0", "#045a8d", "#08306B","#023851",
                                    "grey15", "black"))+
      labs(y = mark, col=input$selectGeneFeature)+
      scale_x_continuous(limits = c(0,600),
                         breaks = c(0, 200,300, 400, 600),
                         labels = c("-3000", "-1000", "TSS", "1000", "3000"))+
      xlab(NULL)
  })
}

