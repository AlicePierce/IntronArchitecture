# -------------------------------------------------------------------------
# File: ui.R
# Pierce et al 2025
# Purpose: Create Scatterplot Shiny App
# -------------------------------------------------------------------------

library(shiny)
library(bslib)

ui <- page_sidebar(
  title="Relationship between Gene Features and Chromatin Features",
  sidebar = sidebar(
    selectInput(
      "selectRegion",
      "Gene Region",
      choices = c("All Gene", "All Exon", "All Intron"),
      selected = "All Gene"
    ),
    selectInput(
      "selectGeneFeature",
      "Gene Feature",
      choices = c("Intron Number", "First Intron Position", "Intron Percent",
                  "First Intron Length", "Gene Length", "Total Exon Length",
                  "Total Intron Length"),
      selected = "Intron Number"
    ),
    selectInput(
      "selectChromatinFeature",
      "Chromatin Feature",
      choices = c("ATAC", "Average Expression", "CG", "CHG", "CHH", "CV",
                  "H2A.W", "H2A.X", "H2A.Z", "H3.1", "H3.3", "H3K14ac",
                  "H3K23ac", "H3K27ac", "H3K27me1", "H3K27me3",
                  "H3K36ac", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3",
                  "H3K56ac", "H3K9ac", "H3K9me1", "H3K9me2", "H4K16ac",
                  "IBM1", "JMJ14", "Pol2", "PolIV", "SDG8", "Tissue Breadth"),
      selected = "ATAC"
    )
  ),
  card(
    card_header(textOutput("dynamicHeader")),
    plotOutput("scatterplot")
  )
)
