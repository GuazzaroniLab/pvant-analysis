This repository contain all data and analysis examples for the manuscript **"Host-Dependent Improvement of GFP Expression in *Pseudomonas putida* KT2440 Using Terminators of Metagenomic Origin"** (now published in [ACS Synthetic Biology]( https://doi.org/10.1021/acssynbio.3c00098)). Details on how the repository is structured and to use the code are provided below.

## Pre-requisites

The original analyses in this work were made using R version 4.1.2 running under [Pop! OS](https://pop.system76.com/) 22.04 LTS. As an attempt to provide reasonable compatibility across R installations, all relevant packages for the analysis are included in an [`renv`](https://rstudio.github.io/renv/articles/renv.html) project environment, but please be sure to be using R version 4.0 or higher when following along. A listing of the packages used in different points during the analysis can be found at the end of the pre-compiled `html` reports found in the examples directory. If you run into any issues running the code described here, please, contact us (see below).

## File structure

The repository is organized as follows:

1.  `data/`: Contains raw data for the flow cytometry and microplate reader assays

    1.1. `flow_cytometry/`: contains directories for the *E. coli* and *P. putida* flow cytometry standard files (.fcs), located in subdirectories for each organism type and experimental batch (named in a DDMMYYYY format).

    Files within these directories can be labeled as `PVANT_BRANCO`/`PSEVA_BRANCO` (corresponding to the empty vectors described in the manuscript), `PVANT SP_BRANCO`/`PSEVA SP_BRANCO` (corresponding to the vectors containing *gfplva* but no promoter) or `PVANT`/`PSEVA_GFP` (corresponding to the vectors containing both *gfplva* and the constitutive promoter).

    1.2. `growth_curves/`: contains reader outputs (.txt) generated by a Victor X3 multilabel microplate reader, and their corresponding plate layouts (.xlsx) for analysis using the ad-hoc R package [mipreadr](github.com/viana-guilherme/mipreadr).

    In the 96-well microplates, in general, different technical replicates are laid in horizontally (column-wise) adjacent wells, while biological replicates are laid in different rows. Samples are labeled with respect to their vector, with `pSEVA`/`pVANT` corresponding to the empty vectors described in the manuscript, `pSEVA.GFP`/`pVANT.GFP` corresponding to the vectors containing *gfplva,* but no promoter, and `pSEVA.P100.GFP`/`pSEVA.P100.GFP` corresponding to the vectors containing both *gfplva* and the constitutive promoter. This information is followed by an underscore and the host (either `#Ecoli` or `#Pputida`), which is also treated as a "blank" in the mipreadr workflow.

2.  `examples/`: Contains step-by-step analysis examples for microplate reader data (microplateAnalysis.Rmd) and flow cytometry data (cytometryAnalysis.Rmd) in RMarkdown format. These documents can also be found as pre-compiled `.html` files that can be read in any browser.

3.  `pvant.RProj`: RStudio project file associated with this working directory and R environment. Please use this file to open RStudio with the correct settings.

4.  `README.md`: Information about this repository (this file).

5.  `utils.R`: Script with custom function definitions (mainly convenience wrappers for handling flow cytometry data and plotting preferences). Used internally in the analysis scripts.

6.  `renv.lock` and `renv/`: Internal files and directory for the package manager

## Running the code

We recommend using RStudio IDE for opening the RMarkdown files available in the examples directory of this repository and following the descriptions provided.

To avoid any potential issues with relative paths for the data directories and to ensure that the correct `renv` used for the analysis is loaded, please, be sure to download the entire repository and, preferably, open the project in RStudio e.g. by double clicking `pvant.Rproj`.

## Contact information

We will be happy to answer if there are any questions or issues found in this repository. Please write to corresponding author María-Eugenia Guazzaroni at [meguazzaroni\@ffclrp.usp.br](mailto:meguazzaroni@ffclrp.usp.br)
