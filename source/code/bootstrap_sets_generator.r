#***************************
#* bootstrap_sets_generator.r
#* Create bootstrap datasets using random resampling with replacement
#* Important notices to keep in mind:
#*    1) stick to folder and data file naming conventions for sets
#*    2) if taxon_A is replaced, R automatically names the first replacement as taxon_A.1 
#*      the second replacement as taxon_A.2 etc. 
#*    3) bootstrap_binary_classifier.r assumes that these conventions are in place.
#*    4) resampling is done separately for each class label "0" and "2", so that bootstrap datasets preserve 
#*       the number of taxa of a given class in the original dataset. 
#*  Cem Ozen, Aug. 9, 2022.  
#**************************

library(geomorph)
library(rgl)
library(ape)
library(phytools)
library(geiger)
library(ggplot2)
library(graphics)
library(NbClust)
library(svgViewR)
library(formattable)
library(htmlwidgets)
library(scatterplot3d)
library(RColorBrewer)
#library(WGCNA)
library(geometry)
library(dispRity)
library(kader)
library( paleotree )
library( HDInterval )
library (MASS)
require( mda )
require( ape )
require( paleotree )
require( nlme )
require( qpcR )
require( klaR )
library(dplyr)
require(sjmisc)
library(data.table)    

runStandAlone       <- FALSE

# INPUT:
baseDir             <- "c:/dino_nature/output/test_bootstrap"                     # base address where individual set folders will be placed 
dataFile            <- "c:/dino_nature/data/femur_compactness_all.csv"            # data file from which bootstrap sets will be generated
boneType            <- "femur"                                                                   # enter "femur" or "ribs" 
applyDataFilter     <- FALSE                                                                   # enter FALSE (no filter) or filter name. (see auxiliary.r) 
Nbootstrapsets      <- 10                                                                        # number of bootstrap sets


create_bootstrap_sets <- function(){
  phyloFDAfile        <- "c:/dino_nature/code/phylo.fda.R"                          # pFDA code
  auxFile             <- "c:/dino_nature/code/auxiliary.r"                          # helper functions
  source(phyloFDAfile)    
  source(auxFile)
  
  setwd(baseDir)
  
  # clear logfiles if they already exist
  to_be_deleted <- list.files(baseDir, pattern = "^logfile")
  if (length(to_be_deleted>0)){
    file.remove(to_be_deleted)
  }
  
  # create log file
  logFile             <- sprintf("logfile_%s.log", Sys.Date())
  cat("Logfile for a session of program: bootstrap_sets_generator.r", sep="\n", file=logFile)
  cat("-----------------------------------------------------", sep='\n',file=logFile, append = T)
  cat("INPUT:", sep="\n", file=logFile, append = T)
  cat(paste0("baseDir: ", baseDir), sep="\n", file=logFile, append = T)
  cat(paste0("boneType: ", boneType), sep="\n", file=logFile, append = T)
  cat("", sep="\n", file=logFile, append = T)
  cat(paste0("Bootstrap data sets have been generated from: ", dataFile), sep="\n", file=logFile, append = T)
  cat(paste0("Applied filters (before set generation) : ", applyDataFilter), sep="\n", file=logFile, append = T)
  
  
  
  data.histo  <- read.csv(dataFile, row.names = 1)                  # main data that will be subject to bootstrapping
  if (applyDataFilter != FALSE){
    data.histo <- manipulateData(data.histo, applyDataFilter)
  }
  
  
  # now separate categories of data (such as diving, non-diving, UNKNOWN").
  # I use the diving.or.not column for this job:
  data0 <- data.histo[data.histo$diving.or.not==0,]   
  data1 <- data.histo[data.histo$diving.or.not==1,]
  data2 <- data.histo[data.histo$diving.or.not==2,]
  data3 <- data.histo[data.histo$diving.or.not=="UNKNOWN",]   # this class will be kept as is, and will be added to every set generated
  
  for (no.set in 1:Nbootstrapsets){
    data0_new <- data0[sample(x=rownames(data0), size = nrow(data0), replace = T),]
    data1_new <- data1[sample(x=rownames(data1), size = nrow(data1), replace = T),]
    data2_new <- data2[sample(x=rownames(data2), size = nrow(data2), replace = T),]
    data.histo.new <- rbind(data0_new, data1_new, data2_new, data3)
    dir.create(sprintf("set_%d", no.set), showWarnings = FALSE) 
    setDataFile <- file.path(baseDir, paste0("set_", no.set), paste0(boneType,"_set_",no.set,".csv"))    # bone data file (with path)
    write.csv(data.histo.new, file = setDataFile)
  }
}

if (runStandAlone){
  result <- create_bootstrap_sets()
}



