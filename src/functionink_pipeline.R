############################
#   functionink_pipeline.R
############################
#
# This script computes the whole functionink pipeline (see
# https://apascualgarcia.github.io/functionInk/). It is assumed
# that the user is following the Vignette described in the documentation and
# hence that the repository was cloned, and that this script is executed with Rstudio.
# The example is set to be executed in its location in the repository and hence
# the paths are automatically set, but you should accomodate these paths to your
# environment.
############################
# USAGE: Provide path and name of the network and of the repository, and set the options desired.
# OUTPUT: A directory located in the same path as the network, containing all output files
############################

#library(here) # you can skip if you use absolute paths
rm(list=ls())
#### START EDITING
# --- Path to repository
pathRepo="10_functionink/functionInk" #  A symbolic link located at 10 pointing to the root of functionink repo 

# --- Name and path of the file containing the network
pathNet="9_predictions/" # absolute path, or one relative to the root of the repo
fileNet="Network_transformation_matrix_Starting2Rep4.tsv"

# ---- STOP EDITING HERE 
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]

# --- Set up the paths and source scripts
pathRepo = paste(this.dir, pathRepo,sep="/")
functionink.scripts=paste(pathRepo,"scripts","analysis_R",sep="/")
setwd(functionink.scripts)
source("extractPartDensity.R")
source("run_pipeline.R")
setwd(this.dir)

# --- Set  path network
pathNet = paste(this.dir, pathNet,sep="/")

# --- Run the pipeline
run_pipeline(fileNet,pathNet = pathNet, pathRepo = pathRepo, # mandatory
weighted=TRUE,directed=FALSE, # options set as in the vignette
types=TRUE,method="Average",mode="internal")
