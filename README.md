# Intron Architecture 

## 🧬 Overview
This repository contains the source code for all analyses and figures associated with the manuscript: *"Intron architecture predicts chromatin features in Arabidopsis thaliana"* 

👉 [Read the preprint here](https://www.biorxiv.org/content/10.1101/2025.10.15.682614v1.full)

## 💽 Data Sources
- Plant Chromatin State Database
- Tissue-Specific Expression
- DNA methylation

## 📁 Repository Structure
```
IntronArchitecture/
├── functions_and_intermediate files/      # R scripts containing analysis and plotting functions
├── figures/                             # Scripts used to generate the manuscript figures
├── shiny_app/                           # Shiny app source code
└── README.md                            # ⭐️ You are here
```

## 🔗 Quick Links
- Functions
- Figure Scripts
- Shiny App

## 💻 Running the Shiny App
To explore the association between intron architecture and chromatin features in Arabidopsis interactively:

1. Download the data from here
2. Download the shiny app folder
3. Install required packages in R
```
# Install required packages
install.packages(c("shiny", "ggplot2", "dplyr", "tidyr", "readr"))
```
4. Run the app locally
```
# Run the app locally
shiny::runApp("shiny_app")
```

## 📝 Citation
Pierce, Alice V., Alan B. Rose, and J. Grey Monroe. 2025. “Intron Architecture Predicts Chromatin Features in Arabidopsis Thaliana.” *bioRxiv*. https://doi.org/10.1101/2025.10.15.682614.

## 🌱 Getting Help
This Repository is maintained by [Alice](https://github.com/AlicePierce) \
For questions or collaborations, feel free to reach out: [Contact Me](mailto:avpierce@ucdavis.edu)

Submit an issue [here](https://github.com/AlicePierce/IntronArchitecture/issues)





