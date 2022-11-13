This repository contain all data and relevant analysis examples for the manuscript **"Termination aware vector design for improved heterologous expression in *Pseudomonas putida* KT2440"** (in preparation). Details on how to use the repository is structured are listed below

## Pre-requisites

The original analyses in this work were made using R version 4.1.2 running under PopOS! 22.04 LTS. All relevant packages used will be included in an renv so there should be reasonable compatibility across operating systems and R installations, but please be sure to be using R version 4.0 or higher when following along. If you run into any issues running the code described here, please, contact us (see below).

## File structure

The repository is organized as follows:

1.  `data/`: Contains raw data for the flow cytometry and microplate reader assays

    1.1. `flow_cytometry/`: contains subdirectories for the E. coli and P. putida experiments, inside which (named in a DDMMYYYY format)

    1.2. `growth_curves/`: contains Victor X3 reading outputs (.txt) and plate maps (.xlsx) for analysis using our ad-hoc R package [mipreadr](github.com/viana-guilherme/mipreadr)

2.  `examples/`: Contains step-by-step example of microplate reader (microplateAnalysis.Rmd) or flow cytometry (cytometryAnalysis.Rmd) analyses in the RMarkdown format.

3.  utils.R: script with custom function definitions (mainly convenience wrappers for handling flow cytometry data and plotting preferences). Used internally in the analysis scripts.

## Running the code

We recommend using RStudio IDE for opening the RMarkdown files available in the examples directory of this repository and following the analysis descriptions provided. To avoid any potential issues with relative paths for the data directories, please, be sure to download the entire repository and open the project in RStudio e.g. by double clicking `pvant.Rproj` or set the working directory accordingly.

## Contact information

We will be happy to answer if there are any questions or issues found in this repository. Please write to María-Eugenia Guazzaroni at [meguazzaroni\@ffclrp.usp.br](mailto:meguazzaroni@ffclrp.usp.br){.email}
