#   loocv_binary.r
#   Leave one out cross validation code suitable for binary classification
#   Cem Ozen, Aug. 10, 2022

#   STACKHEAD 8 
#
#   Run a loop over every member of the training set. 
#   For each species in the training set you do the following
#       i. Remove the species from the training set.
#      ii. Run the analysis, with the removed species as the test species.
#     iii. This should give you the set of 100 predictions – and we look at the median to get the net prediction
#   The LOOCV tells whether their method can correctly predict the animals in the training set.   

library(parallel)
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

runStandAlone       <- FALSE     # set TRUE if this file is used as main code, otherwise FALSE for using it as a module

# provide user input
baseDir             <- "c:/dino_nature/output/test_loocv"                         # output directory
boneType            <- "femur"                                                                    # enter "femur" or "ribs"                             
dataFile            <- "c:/dino_nature/data/femur_compactness_all.csv"        # bone data file (with path)
applyDataFilter     <- "F0D02"                                                                      # enter FALSE (no filter) or filter name. (see auxiliary.r) 
addTreeTips         <- FALSE                                                                        # FALSE for LOOCV normally (unless still including test taxa in which new variants have been added)
Ntrees              <- 100                                                                          # no of random trees to be generated
calculateTrees      <- TRUE                                                                         # calculate (TRUE) or read (FALSE) trees
saveTrees           <- FALSE                                                                        # save (TRUE) trees (if they are calculated here)



MC_loocv <- function(set_no, logFile, baseDir, auxFile, phyloFDAfile, Ntrees, data.histo.keep, data.histo.test, scaled.phylos){
  
  message <- paste0("******* Working on LOOCV set:", set_no, "********")
  print(message)
  cat(message, sep="\n", file=logFile, append = T)
  
  source(phyloFDAfile)    
  source(auxFile) 
  
  
   setwd(baseDir)
   dir.create(paste0("set_", set_no), showWarnings = FALSE)
   workDir <- file.path(baseDir, paste0("set_", set_no))
   setwd(workDir)
   
  # # reset current loocv data to our kept data consisting of only training taxa with their genuine labels
   data.histo <- data.histo.keep
   # set current taxon's diving.or.not to UNKNOWN
   data.histo[set_no, "diving.or.not"] <- "UNKNOWN"
  
   # we can just make a prediction for the current taxon in the loocv. But we can also still predict the real test taxa that
   # separated out above. we will add them to the data now to predict them as well (OPTIONAL)
   data.histo <- rbind(data.histo, data.histo.test)    # only taxa w/ UNKNOWN value for diving.or.not in data are test taxa
   
   
   # our test set is made of the current taxon
   test.set <- rownames(data.histo[data.histo$diving.or.not == "UNKNOWN", ])
   test.set
   
   
   # crucial step: order taxa in data after tree order. first ensure data and tree are in concord
   # (they are made of the exact same list of taxa)
   random.num <- round(runif(1, min = 1, max = Ntrees), digits = 0)
   tree.temp <- scaled.phylos[[random.num]]
   concord <- concord_tree_data(tree.temp, data.histo)
   tree.temp <- concord$tree
   data.histo <- concord$data
   data.histo <- data.histo[tree.temp$tip.label,]
   
   
   # prepare the input vectors and label vectors for classification
   diving.or.not.all <-data.histo[, "diving.or.not"]; names(diving.or.not.all) <- row.names(data.histo)
   compactness <- data.histo[ , "Global.compactness" ]; names(compactness) <- row.names(data.histo)
   diameter <- log10(data.histo[ , "MD..mm." ]); names(diameter) <- row.names(data.histo)
   
   
   variable = diving.or.not.all
   test.taxa <- which(variable == "UNKNOWN", useNames = TRUE)
   training.set <- names(variable [-test.taxa])
   test.set <- names(variable[test.taxa])
   names(test.taxa) <- test.set
   training.compactness <- compactness [training.set]
   training.diameter <- diameter [training.set]
   training.groups   <- diving.or.not.all [training.set]
   XA <- as.matrix(cbind(compactness, diameter))
   X <- as.matrix(cbind(training.compactness, training.diameter))
   
   # tree.temp which is in concord with data (combination of training and test)
   # but we also need tree.training from that by pruning test set taxa
   drop.taxa <- tree.temp$tip.label[ ! tree.temp$tip.label %in% training.set]
   tree.training <- drop.tip( tree.temp, drop.taxa)
   
   filename_stem <- "ecology"
   opt1<-optLambda(X,training.groups,tree.training,idc=filename_stem )
   lambda <- as.numeric(opt1$optlambda) # Replace with the optimal lambda value from above.
   #cat(paste0("* lambda (calculated from a random tree) = : ", lambda), sep="\n", file=logFile, append = T)
   print(lambda)
   
   # important tests: to prevent mismatching between tree tips and data or derived objects, data must have been ordered after tree so far.
   # this should be separately true for combined data (training + test) and its tree and true for training data and its tree.
   # this test ensures that we can address any taxa either through its index or its name, without any conflict.
   test1 <- all(names(compactness) == tree.temp$tip.label) & all(names(diameter) == tree.temp$tip.label) & all(rownames(XA) == tree.temp$tip.label) & all(names(variable) == tree.temp$tip.label)
   test2 <- all(names(training.compactness) == tree.training$tip.label) & all(names(training.diameter) == tree.training$tip.label) & all(rownames(X) == tree.training$tip.label) & all(names(training.groups) == tree.training$tip.label)
   test3 <- all(tree.temp$tip.label[as.vector(test.taxa)] == names(test.taxa))
   if (!(test1 & test2 & test3)){
     message <- "* ERROR: check index-name alignment!"
     cat(message, sep="\n", file=logFile, append = T)
     print(message)
     stop("* ERROR: Conflicts exist between some objects. Check integer indices and names of taxa!")
   } else {
     message <- "* index-name alignment check: OK"
     cat(message, sep="\n", file=logFile, append = T)
     print(message)
   }

  pfdas<-list()
  test.lambda.results <- list()
  predicted.ecologies <- matrix(data = NA, nrow = length(test.taxa), ncol = Ntrees)
  rownames(predicted.ecologies) <- names(test.taxa)
  rate.accuracy.subaequous <- c()
  rate.accuracy.non.subaequous<- c()
  probability.subaequous <-matrix(data = NA, nrow = length(test.taxa), ncol = Ntrees)
  rownames(probability.subaequous) <- names(test.taxa)
  probability.non.subaequous <- matrix(data = NA, nrow = length(test.taxa), ncol = Ntrees)
  rownames(probability.non.subaequous) <- names(test.taxa)

  boundaryline_coeffs <- matrix(data = NA, nrow = Ntrees , ncol = 3)   # storage for boundary line coefficients
  colnames(boundaryline_coeffs) <- c("coeff_diameter","coeff_compactness","coeff_intercept")

  confusion.matrix <- matrix(data = NA, nrow = Ntrees, ncol = 4)  # to store confusion matrix of training for all trees

  for( i in 1:length( scaled.phylos ) ) {
    tree.temp <- scaled.phylos[[i]]
    # from data.histo we had already removed taxa that are not in a randomly sampled tree above. we only need to remove taxa from each random
    # tree. we do not need worry about taxa order as data.histo was already ordered in a way to match the order of the tree after its tips
    # are rid of unwanted taxa:
    drop.taxa <- tree.temp$tip.label[ ! tree.temp$tip.label %in% rownames(data.histo)]
    tree.temp <- drop.tip(tree.temp , drop.taxa)
    pfda <- phylo.fda.pred(XA,variable,tree.temp$tip.label,tree.temp,test.taxa,val = lambda, eqprior=TRUE, priin = 1)
    pfdas [[i]]<-pfda

    boundaryline_coeffs[[i,"coeff_intercept"]]  <- pfda$fit$coefficients[[1]]
    boundaryline_coeffs[[i,"coeff_compactness"]]  <- pfda$fit$coefficients[[2]]
    boundaryline_coeffs[[i,"coeff_diameter"]]  <- pfda$fit$coefficients[[3]]

    #compile rate of missclasificacion
    true.non.diving <- pfda$confusion[1,1]
    total.non.diving <- pfda$confusion[1,1] + pfda$confusion[2,1]
    rate.accuracy.non.diving <- true.non.diving / total.non.diving
    rate.accuracy.non.subaequous[i]<-rate.accuracy.non.diving

    true.diving <- pfda$confusion[2,2]
    total.diving <- pfda$confusion[2,2] + pfda$confusion[1,2]
    rate.accuracy.diving <- true.diving / total.diving
    rate.accuracy.subaequous[i] <- rate.accuracy.diving

    test.class <- as.matrix(predict(pfda, pfda$DATAtest, type="class"))
    predicted.ecologies[,i]<- test.class[,1]

    # Discriminant scores scores for the fossils are retrieved this way.
    test.variates <- round(predict(pfda, pfda$DATAtest, type="variates"), digits = 2)

    # Posterior probabilities of group affiliations are available, too.
    test.prob <- round(predict(pfda, pfda$DATAtest, type="posterior"), digits = 2)
    probability.non.subaequous[, i] <- test.prob[,1]
    probability.subaequous[, i] <- test.prob[,2]

    # And now let's put everything together (here shown for the fossils):
    test.results <- cbind(test.class, test.prob, test.variates)
    colnames(test.results) <- c("predicted class", "P(not-diving)", "P(diving)", "DA1")
    rownames(test.results) <- names(test.taxa)

    test.lambda.results[[i]]<- test.results

    # confusion matrix
    confusion.matrix[[i,1]] <- pfda$confusion[1,1]
    confusion.matrix[[i,2]] <- pfda$confusion[1,2]
    confusion.matrix[[i,3]] <- pfda$confusion[2,1]
    confusion.matrix[[i,4]] <- pfda$confusion[2,2]
  }

  write.csv(boundaryline_coeffs, file = "./BOUNDARY_LINE_COEFFS.csv")

  #number of times predicted as subaequous forager for the 15 taxa that we are to predict
  predictions<-row_count(as.data.frame(predicted.ecologies), count = "2")
  predictions.final<- predictions$rowcount ; names(predictions.final) <- names(test.taxa)
  predictions.final


  #get rates accuracy by ecological category
  rates.accuracy <- cbind(median(rate.accuracy.non.subaequous), median(rate.accuracy.subaequous)) ; colnames(rates.accuracy) <- c("Rate accuracy terrestrial","Rate accuracy subaequous forager")
  rates.accuracy

  #summarise median and maybe range values for predictions of both categories
  medians.non.subaequous <- rowMedians(probability.non.subaequous)
  medians.subaequous <- rowMedians(probability.subaequous)
  medians.subaequous          # median probability of diving (over the given number of trees)


  variable[test.set]           # original labels of the test set
  probability.subaequous       # probabilities of diving from each tree
  write.csv(probability.subaequous, file = "./PROBABILITIES.csv")


  predictions.final            # number of times test taxa are predicted to be diving (over the given number of trees)
  summary.results <- matrix(ncol=4, nrow=length(test.set))
  for( i in 1:length(test.set) ) {
    all.median.probabilities <- c(medians.non.subaequous[[test.set[[i]]]], medians.subaequous[[test.set[[i]]]])
    predicted.class.index <- which.max( all.median.probabilities )     # which returns the index of max
    predicted.class <- colnames(pfda$confusion)[[predicted.class.index]]
    summary.results[[i, 1]] <- data.histo.keep[test.set[[i]],"diving.or.not"]    # test taxa original labels
    summary.results[[i, 2]] <- predicted.class                     # test taxa predicted labels
    summary.results[[i, 3]] <- predictions.final[[test.set[[i]]]]  # number of times test taxa are predicted to be diving (over the given number of trees)
    summary.results[[i, 4]] <- medians.subaequous[[test.set[[i]]]] # median probability of diving (over the given number of trees)
  }
  colnames(confusion.matrix) <- c("[1,1]", "[1,2]", "[2,1]", "[2,2]")
  write.csv(confusion.matrix, file = "./CONFMAT.csv")
  confusion.matrix
  colnames(summary.results) <- c("actual class", "predicted class", "N", "median P(diving)")
  rownames(summary.results) <- test.set
  write.csv(summary.results, file = "./SUMMARY.csv")
  summary.results
  #cat("", sep="\n", file=logFile, append = T)
}

run_loocv <- function(){
  workDir <- baseDir   
  setwd(workDir)
  
  phyloFDAfile        <- "c:/dino_nature/code/phylo.fda.R"                          # pFDA code
  auxFile             <- "c:/dino_nature/code/auxiliary.r"                          # helper functions
  source(phyloFDAfile)    
  source(auxFile) 
  
  # clear logfiles if they already exist
  to_be_deleted <- list.files(workDir, pattern = "^logfile")
  if (length(to_be_deleted>0)){
    file.remove(to_be_deleted)
  }
  # create log file
  logFile             <- file.path(baseDir, sprintf("logfile_%s.log", Sys.Date()))
  cat("Logfile for a session of program: loocv_binary.r", sep="\n", file=logFile)
  cat("-----------------------------------------------------", sep='\n',file=logFile, append = T)
  cat("", sep="\n", file=logFile, append = T)
  cat("WARNING: When comparing estimated labels to true labels, be aware of how the data was treated", sep="\n", file=logFile, append = T)
  cat("         (the label column, any filters applied, etc.)", sep="\n", file=logFile, append = T)
  cat("", sep="\n", file=logFile, append = T)
  cat("INPUT:", sep="\n", file=logFile, append = T)
  cat(paste0("baseDir: ", baseDir), sep="\n", file=logFile, append = T)
  cat(paste0("phyloFDAfile: ", phyloFDAfile), sep="\n", file=logFile, append = T)
  cat(paste0("auxFile: ", auxFile), sep="\n", file=logFile, append = T)
  cat(paste0("boneType: ", boneType), sep="\n", file=logFile, append = T)
  cat(paste0("dataFile: ", dataFile), sep="\n", file=logFile, append = T)
  cat(paste0("applyDataFilter: ", applyDataFilter), sep="\n", file=logFile, append = T)
  cat(paste0("addTreeTips: ", addTreeTips), sep="\n", file=logFile, append = T)
  cat(paste0("Ntrees: ", Ntrees), sep="\n", file=logFile, append = T)
  cat(paste0("calculateTrees: ", calculateTrees), sep="\n", file=logFile, append = T)
  cat(paste0("saveTrees: ", saveTrees), sep="\n", file=logFile, append = T)
  cat("", sep="\n", file=logFile, append = T)
  cat("", sep="\n", file=logFile, append = T)
  cat("OPERATIONAL CHECKS:", sep="\n", file=logFile, append = T)
  
  
  # duplicated taxa: if present in data but not tree, we need to get the (original taxa, duplicate) pairs 
  if (addTreeTips){
    data.histo <- read.csv(dataFile, row.names = 1)
    base_duplicate_pairs <- get_base_duplicate_pairs(data.histo)
    message <- "* Regression data contains multiple copies of some taxa (user added variants of existing test taxa \n   or replacements due to bootstrap are present) which are not in the phylo tree:"
    print(message)
    cat(message, sep="\n", file=logFile, append = T)
    cat("     original taxon,   duplicate   ", sep="\n", file=logFile, append = T)
    write.table(base_duplicate_pairs, file=logFile, append = T, sep = ",", row.names = T, col.names = F)
    cat(" ", sep="\n", file=logFile, append = T)  
    print(base_duplicate_pairs)
  }
  
  # read phylogenetic tree and supplementary data to manipulate the tree (drop unwanted taxa, calibrate nodes, ultrametricize) then generate
  # random samples of scaled trees
  if (boneType == "femur"){
    data.histo     <- read.csv("c:/dino_nature/data/femur_compactness_all.csv", row.names = 1)
    original.phylo <- read.nexus("c:/dino_nature/data/tree_femur_final.nex")
    strata.times   <- read.csv("c:/dino_nature/data/femur_strata_multiphylo.csv",  header = T, row.names = 1)
    taxon.strata   <- read.csv("c:/dino_nature/data/femur_taxa_multiphylo.csv",  header = T, row.names = 1) 
    calibrated.nodes <-	read.csv("c:/dino_nature/data/tree_femur_nodes_final.csv", sep = ",", dec = ".", stringsAsFactors = F,colClasses = c("NODE.NUMBER" ="numeric", "MINIMUM.AGE"="numeric"), header = T, row.names = 1)
  } else {      # for ribs
    data.histo     <- read.csv("c:/dino_nature/data/ribs_compactness_all.csv", row.names = 1)
    original.phylo <- read.nexus("c:/dino_nature/data/tree_ribs_final.nex")
    strata.times   <- read.csv("c:/dino_nature/data/ribs_strata_multiphylo.csv",  header = T, row.names = 1)
    taxon.strata   <- read.csv("c:/dino_nature/data/ribs_taxa_multiphylo.csv",  header = T, row.names = 1) 
    calibrated.nodes <-	read.csv("c:/dino_nature/data/tree_ribs_nodes_final.csv", sep = ",", dec = ".", stringsAsFactors = F,colClasses = c("NODE.NUMBER" ="numeric", "MINIMUM.AGE"="numeric"), header = T, row.names = 1)
  }
  
  # before random scaled trees are generated we need to add tips to our main tree for all duplicated taxa present in regression data.
  if (addTreeTips){
    # we can only add a tip for a duplicated taxon if the original copy of that taxon is present in the tree. It is possible that some taxa may be present in the 
    # data that do not exist in the main tree. In bootstrap dataset generation, it is possible that such a taxon is duplicated, and may therefore be present in
    # base_duplicate_pairs. We need to eliminate such rows for they would cause and error as we cannot add a tip for them. 
    mask <- base_duplicate_pairs[,"bases"] %in% original.phylo$tip.label
    base_duplicate_pairs <- base_duplicate_pairs[mask,]
    # for each duplicate add a tip at the node its base is attached
    for (i in 1:nrow(base_duplicate_pairs)){
      original.phylo <- phyloDuplicateTip(original.phylo, base_duplicate_pairs[[i,1]], base_duplicate_pairs[[i,2]])
    }
    # check which taxa the added tips share the final node with
    for (i in 1:nrow(base_duplicate_pairs)){
      print(findSpeciesWithCommonRoot(original.phylo,base_duplicate_pairs[[i,1]]))
    }
    # now add duplicate entries to taxon.strata
    for (i in 1:nrow(base_duplicate_pairs)){
      taxon.strata[base_duplicate_pairs[[i,2]],] <- taxon.strata[base_duplicate_pairs[[i,1]],]
    }
    # now add duplicate entries to data needed for tree operations
    for (i in 1:nrow(base_duplicate_pairs)){
      data.histo[base_duplicate_pairs[[i,2]],] <- data.histo[base_duplicate_pairs[[i,1]],]
    }
    message <- "* Phylo tree and the main data and strata files have been modified to include duplicated taxa before the scaled random trees are generated."
    print(message)
    cat(message, sep="\n", file=logFile, append = T)
  }
  
  calibrated.nodes <- calibrated.nodes[, 2]
  taxon.strata.list<- list(strata.times, taxon.strata)
  
  # check ultrametricity
  if (!is.ultrametric(original.phylo)){
    original.phylo <- force.ultrametric(original.phylo)
  }   
  # drop from tree all taxa in tree but not in data 
  drop.taxa <- original.phylo$tip.label[ ! original.phylo$tip.label %in% rownames(data.histo)]
  original.phylo <- drop.tip(original.phylo, drop.taxa)
  if (calculateTrees){
    scaled.phylos <-bin_timePaleoPhy (original.phylo, timeList = taxon.strata.list,
                                      type = "mbl", vartime = 1, ntrees = Ntrees,
                                      nonstoch.bin = FALSE, randres = T, 
                                      timeres = F, sites = NULL, 
                                      point.occur = FALSE, add.term = T, 
                                      inc.term.adj = T, dateTreatment = "firstLast", 
                                      node.mins = calibrated.nodes, noisyDrop = TRUE, plot = F)
    if (saveTrees){    # save trees if need be
      for( i in 1:Ntrees ){
        write.nexus(scaled.phylos[[i]], file = sprintf("./tree_%d.nex", i))
      }
    } 
  } else {    # trees are not calculated, then read them in
    for( i in 1:Ntrees ){
      scaled.phylos[[i]] <- read.nexus(sprintf("./tree_%d.nex", i))
    }
  }
  
  
  # regression data
  data.histo <- read.csv(dataFile, row.names = 1) 
  # apply filters, if any (selecting taxa based on Flying or Diving labels or some other column)
  if (applyDataFilter != FALSE){
    data.histo <- manipulateData(data.histo, applyDataFilter)
  }
  
  # we will first separate out the real test taxa (OPTIONAL)
  test.taxa.patterns <- c("^Spinosaurus", "^Baryonyx", "^Suchomimus")
  test.set.indices <- lapply(test.taxa.patterns, function(x) grep(x, rownames(data.histo)))
  test.set.indices <- unlist(test.set.indices)
  test.set <- rownames(data.histo[test.set.indices,])
  data.histo.test <- data.histo[test.set,] 
  
  
  
  # In LOOCV analysis, one normally takes the training set (Nt taxa) and using Nt-1 of the set as training, one make a prediction for the label of the 
  # left out member and repeat this over the whole set. The tarining part of the data is:
  data.histo <- data.histo[data.histo$diving.or.not != "UNKNOWN",]
  # drop any taxon that is not in main tree, in order not to make such a taxon an LOOCV test case which would create errors since they are not in the tree
  mask <- rownames(data.histo) %in% original.phylo$tip.label
  data.histo <- data.histo[mask,]
  
  message <- "* Selected categorical labels column: diving.or.not"
  cat(message, sep="\n", file=logFile, append = T)
  print(message)
  # important check about the column we will use as array of categorical labels 
  test <- all(data.histo$Diving == data.histo$diving.or.not) 
  if (!test) {
    message <- "* WARNING: diving.or.not column is NOT identical to Diving column!"
    cat(message, sep="\n", file=logFile, append = T)
    print(message)
  } else {
    message <- "* diving.or.equal column is identical to Diving column."
    cat(message, sep="\n", file=logFile, append = T)
    print(message)
  }
  test <- all(data.histo$Flying == 0) # check if data contains flying taxa 
  if (!test) {
    message <- "* Regression includes flying taxa."
    cat(message, sep="\n", file=logFile, append = T)
    print(message)
  } else {
    message <- "* Regression excludes flying taxa."
    cat(message, sep="\n", file=logFile, append = T)
    print(message)
  }
  # check categorical labels to be considered in regression
  diving.classes <- levels(factor(data.histo$diving.or.not)) 
  message <- sprintf("* Diving labels that are present in regression: %s", paste(diving.classes, collapse = ','))
  cat(message, sep="\n", file=logFile, append = T)
  print(message)
  if (length(diving.classes) != 2){
    message <- "* ERROR: this code works with 2 categorical labels but that's not the case with this data. Check the data!"
    cat(message, sep="\n", file=logFile, append = T)
    stop(message)
  }
  
  # now let's keep our data (which contains only training part with their categorical labels intact). In the LOOCV loop we will switch the labels to UNKNOWN
  # one taxon at a time and estimate its label as it is normally done in binary_classifier.r
  data.histo.keep <- data.histo
  
  N.loocv.sets <- nrow(data.histo.keep)
  
  #------------- Monte Carlo runs in series -------------------
  
  #system.time(loocv_output <- lapply(1:N.loocv.sets, MC_loocv, baseDir=baseDir, auxFile=auxFile, phyloFDAfile=phyloFDAfile, Ntrees=Ntrees, logFile=logFile, data.histo.keep=data.histo.keep, scaled.phylos=scaled.phylos))
  
  #------------- Monte Carlo runs in parallel ------------------
  
  numCores <- detectCores()   # how many cores do we have?
  cl <- makeCluster(numCores) # launch a cluster of processes
  
  clusterEvalQ(cl, library(geomorph))
  clusterEvalQ(cl, library(lme4))
  clusterEvalQ(cl, library(geomorph))
  clusterEvalQ(cl, library(rgl))
  clusterEvalQ(cl, library(ape))
  clusterEvalQ(cl, library(phytools))
  clusterEvalQ(cl, library(geiger))
  clusterEvalQ(cl, library(ggplot2))
  clusterEvalQ(cl, library(graphics))
  clusterEvalQ(cl, library(NbClust))
  clusterEvalQ(cl, library(svgViewR))
  clusterEvalQ(cl, library(formattable))
  clusterEvalQ(cl, library(htmlwidgets))
  clusterEvalQ(cl, library(scatterplot3d))
  clusterEvalQ(cl, library(RColorBrewer))
  #clusterEvalQ(cl, library(WGCNA))
  clusterEvalQ(cl, library(geometry))
  clusterEvalQ(cl, library(dispRity))
  clusterEvalQ(cl, library(kader))
  clusterEvalQ(cl, library( paleotree ))
  clusterEvalQ(cl, library( HDInterval ))
  clusterEvalQ(cl, library (MASS))
  clusterEvalQ(cl, require( mda ))
  clusterEvalQ(cl, require( ape ))
  clusterEvalQ(cl, require( nlme ))
  clusterEvalQ(cl, require( qpcR ))
  clusterEvalQ(cl, require( klaR ))
  clusterEvalQ(cl, library(dplyr))
  clusterEvalQ(cl, require(sjmisc))
  clusterEvalQ(cl, library(data.table))
  
  
  loocv_output <- parLapply(cl, 1:N.loocv.sets, MC_loocv, baseDir=baseDir, auxFile=auxFile, phyloFDAfile=phyloFDAfile, Ntrees=Ntrees, logFile=logFile, data.histo.keep=data.histo.keep, data.histo.test=data.histo.test, scaled.phylos=scaled.phylos)
  stopCluster(cl)
}

if (runStandAlone){
  result <- system.time(run_loocv())
  write.table(as.matrix(result), file=file.path(baseDir, "runtime.log"), row.names=T, col.names=F)
}







