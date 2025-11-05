# Intron Architecture

## 🧬 Overview

This repository contains the source code for all analyses and figures associated with the manuscript: `"Intron architecture predicts chromatin features in Arabidopsis thaliana"`

👉 [Read the preprint here](https://www.biorxiv.org/content/10.1101/2025.10.15.682614v1.full)

## 🔗 Quick Links

-   [Shiny App: Relationship between intron architecture and chromatin features](https://github.com/AlicePierce/IntronArchitecture/scatterplots_shiny)
-   [Shiny App: Distribution of chromatin features binned by intron features](https://github.com/AlicePierce/IntronArchitecture/metaplots_shiny)
-   [Instructions for downloading and running the Shiny Apps](https://github.com/AlicePierce/IntronArchitecture?tab=readme-ov-file#-running-the-shiny-app)
-   [Functions](https://github.com/AlicePierce/IntronArchitecture/functions_and_intermediates)
-   [Figure Scripts](https://github.com/AlicePierce/IntronArchitecture/figures)

## ⚙️ Installation

NOTE: Installing/Cloning this entire repository is not required to run the Shiny apps. To only run the apps, follow the provided [instructions for downloading and launching them](https://github.com/AlicePierce/IntronArchitecture?tab=readme-ov-file#-running-the-shiny-app).

To run the code maintained in this repository:

1.  Clone this repository:

```         
git clone https://github.com/AlicePierce/IntronArchitecture.git
cd IntronArchitecture
```

2.  Install the following dependencies:\
    This repository was developed in R (≥ 4.2.0) and uses several packages from CRAN and Bioconductor. Use the following code to automatically install all required packages:

``` r
# List of required packages
packages <- c(
  "data.table", "tidyverse", "rtracklayer", "Hmisc", "reshape2",
  "randomForest", "Biostrings", "GenomicRanges", "RColorBrewer",
  "Rsamtools", "writexl", "cowplot", "shiny", "bslib"
)

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install missing packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch(
      install.packages(pkg, dependencies = TRUE),
      error = function(e) BiocManager::install(pkg, ask = FALSE)
    )
  }
}
```

3.  Download the following files and place them in the `Input-files` folder before running the code:

-   TAIR10 Genome Release: [TAIR10_GFF3_Genes.gff](https://www.arabidopsis.org/download/list?dir=Genes%2FTAIR10_genome_release%2FTAIR10_gff3)
-   TAIR10 FASTA: [TAIR10_chr_all.fas.gz](https://www.arabidopsis.org/download/list?dir=Genes%2FTAIR10_genome_release%2FTAIR10_chromosome_files)
-   Arabidopsis Representative Gene Models: [TAIR10_representative_gene_models](https://www.arabidopsis.org/download/list?dir=Genes%2FTAIR10_genome_release%2FTAIR10_gene_lists)
-   Chromatin Features (BW format) from [Yue Liu et al. 2018](https://academic.oup.com/nar/article/46/D1/D1157/4429024)
-   Cytosine Methylation - follow instructions from [Monroe et al. 2022](https://www.nature.com/articles/s41586-021-04269-6)
-   Tissue-specific Expression from [Mergner et al. 2020](https://www.nature.com/articles/s41586-020-2094-2)
-   Gene Duplicates from [Kenchanmane Raju et al. 2023](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.19227)
-   Intron Coordinates for Gene Duplicate Analysis from [Monroe et al. 2022](https://www.nature.com/articles/s41586-021-04269-6) : [41586_2021_4269_MOESM4_ESM.txt](https://www.nature.com/articles/s41586-021-04269-6#Sec45) (supplementary data 2)

4.  Download the Shiny App Inputs from [here](https://doi.org/10.5281/zenodo.17536897) and place the files in their respective Shiny App folders
5.  Create `Input-files/` and `Output-files/` folders

## 📁 Repository Structure

```         
IntronArchitecture/
├── functions_and_intermediates/         # R scripts containing analysis and plotting functions
├── figures/                             # Scripts used to generate the manuscript figures
├── scatterplots_shiny/                  # Shiny app source code
├── metaplots_shiny/                     # Shiny app source code
├── Input-files/                         # Input files
├── Output-files/                        # Output files
├── LICENSE                              # MIT License
└── README.md                            # ⭐️ You are here
```

## 💻 Running the Shiny Apps

### Dependencies

Install the required packages:

``` r
# Install required packages
install.packages(c("shiny", "bslib", "ggplot2"))
```

Download the required data from Zenodo DOI: [10.5281/zenodo.17536896](https://doi.org/10.5281/zenodo.17536896)

Upload the data to the corresponding Shiny App folder and Run the App

### Scatter Plots Shiny App

To explore the relationship between gene features and chromatin features interactively:

1.  Download the [shiny app folder](https://github.com/AlicePierce/IntronArchitecture/scatterplots_shiny)
2.  Download [ScatterData.rds](https://zenodo.org/records/17536897) and place it in the `scatterplots_shiny` folder
3.  Run the app locally using the following code:

``` r
# Run the app locally
shiny::runApp("scatterplots_shiny")
```

### Metaplots Shiny App

To explore the distribution of chromatin features binned by gene architecture interactively:

1.  Download the [shiny app folder](https://github.com/AlicePierce/IntronArchitecture/metaplots_shiny)
2.  Download [metaplotsData.rds](https://zenodo.org/records/17536897) and place it in the `metaplots_shiny` folder
3.  Run the app locally using the following code:

``` r
# Run the app locally
shiny::runApp("metaplots_shiny")
```

## 📝 Citation

Alice V. Pierce, Alan B. Rose, and J. Grey Monroe. 2025. "Intron Architecture Predicts Chromatin Features in Arabidopsis Thaliana." *bioRxiv*. <https://doi.org/10.1101/2025.10.15.682614>.

## 🌱 Getting Help

This repository is maintained by [Alice](https://github.com/AlicePierce)\
For questions or collaborations, feel free to reach out: ✉️ [Contact Me](mailto:avpierce@ucdavis.edu)\
Submit an issue [here](https://github.com/AlicePierce/IntronArchitecture/issues)
