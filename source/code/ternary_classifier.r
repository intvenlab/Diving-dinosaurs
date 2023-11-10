# binary classifier with tree modifications. 
# Cem Ozen, Aug 8. 2022 

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


runStandAlone      <- FALSE

# provide user input
baseDir             <- "c:/dino_nature/output/D9_10_27_2022/single"                      # output directory
boneType            <- "femur"                                                                # enter "femur" or "ribs" 
dataFile            <- "c:/dino_nature/data/D9_femur.csv"         # bone data file (with path)
applyDataFilter     <- FALSE                                                                   # enter FALSE (no filter) or filter name. (see auxiliary.r) 
addTreeTips         <- TRUE                                                                      # if taxa duplicates are present in reg. data but not tree, TRUE, otherwise FALSE
Ntrees              <- 100                                                                       # no of random trees to be generated
calculateTrees      <- TRUE                                                                      # calculate (TRUE) or read (FALSE) trees
saveTrees           <- FALSE                                                                     # save (TRUE) trees (if they are calculated here)


run_single <- function(){
  setwd(baseDir)
  phyloFDAfile        <- "c:/dino_nature/code/phylo.fda.R"                          # pFDA code
  auxFile             <- "c:/dino_nature/code/auxiliary.r"                          # helper functions
  source(phyloFDAfile)    
  source(auxFile) 
  
  # clear logfiles if they already exist
  to_be_deleted <- list.files(baseDir, pattern = "^logfile")
  if (length(to_be_deleted>0)){
    file.remove(to_be_deleted)
  }
  # create log file
  logFile             <- sprintf("logfile_%s.log", Sys.Date())
  cat("Logfile for a session of program: binary_classifier.r", sep="\n", file=logFile)
  cat("-----------------------------------------------------", sep='\n',file=logFile, append = T)
  cat("INPUT:", sep="\n", file=logFile, append = T)
  cat(paste0("baseDir: ", baseDir), sep="\n", file=logFile, append = T)
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
  
  #message <- paste0("* Number of taxa (after application of filters, if any) = ",nrow(data.histo))
  #cat(message, sep="\n", file=logFile, append = T)
  #print(message)
  
  
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
  # eliminate "UNKNOWN", it is not a categorical label but indicator for unknown category
  diving.classes <- diving.classes[!diving.classes %in% c("UNKNOWN")]
  message <- sprintf("* Diving labels that are present in regression: %s", paste(diving.classes, collapse = ','))
  cat(message, sep="\n", file=logFile, append = T)
  print(message)
  if (length(diving.classes) != 3){
    message <- "* ERROR: this code works with 3 categorical labels but that's not the case with this data. Check the data!"
    cat(message, sep="\n", file=logFile, append = T)
    stop(message)
  }
  
  
  test.taxa.patterns <- c("^Spinosaurus", "^Baryonyx", "^Suchomimus")
  test.set.indices <- lapply(test.taxa.patterns, function(x) grep(x, rownames(data.histo)))
  test.set.indices <- unlist(test.set.indices)
  test.set <- rownames(data.histo[test.set.indices,])
  test.set
  
  
  # remove taxa with UNKNOWN labels if they are not in our test set
  data.histo.test <- data.histo[test.set,] 
  unknown <-which(data.histo$diving.or.not == "UNKNOWN", useNames = TRUE)
  data.histo <- data.histo[-unknown,]
  data.histo <- rbind(data.histo, data.histo.test)    # only taxa w/ UNKNOWN value for diving.or.not in data are test taxa
  data.histo[data.histo$diving.or.not == "UNKNOWN", ] # double-check that only UNKNOWN labels are the test taxa. 
  
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
  
  #message <- paste0("* Number of training taxa  = ",length(tree.training$tip.label))
  #cat(message, sep="\n", file=logFile, append = T)
  #print(message)
  
  filename_stem <- "ecology" 
  opt1<-optLambda(X,training.groups,tree.training,idc=filename_stem )
  lambda <- as.numeric(opt1$optlambda) # Replace with the optimal lambda value from above.
  cat(paste0("* lambda (calculated from a random tree) = : ", lambda), sep="\n", file=logFile, append = T)
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
  rate.accuracy.terrestrial <- c()
  rate.accuracy.amphibious <- c()
  rate.accuracy.aquatic <- c()
  rate.accuracy.overall <- c()
  
  probability.terrestrial <- matrix(data = NA, nrow = length(test.taxa), ncol = Ntrees)   # 0 <> terrestrial
  rownames(probability.terrestrial) <- names(test.taxa) 
  probability.amphibious <-matrix(data = NA, nrow = length(test.taxa), ncol = Ntrees)     # 1 <> amphibious
  rownames(probability.amphibious) <- names(test.taxa) 
  probability.aquatic <-matrix(data = NA, nrow = length(test.taxa), ncol = Ntrees)        # 2 <> aquatic
  rownames(probability.aquatic) <- names(test.taxa) 
  
  boundaryline1_coeffs <- matrix(data = NA, nrow = Ntrees , ncol = 3)   # storage for boundary line coefficients
  colnames(boundaryline1_coeffs) <- c("coeff_diameter","coeff_compactness","coeff_intercept")
  boundaryline2_coeffs <- matrix(data = NA, nrow = Ntrees , ncol = 3)   # storage for boundary line coefficients
  colnames(boundaryline2_coeffs) <- c("coeff_diameter","coeff_compactness","coeff_intercept")
  
  confusion.matrix <- matrix(data = NA, nrow = Ntrees, ncol = 9)  # to store confusion matrix of training for all trees 

  for( i in 1:length( scaled.phylos ) ) { 
    tree.temp <- scaled.phylos[[i]]    
    # from data.histo we had already removed taxa that are not in a randomly sampled tree above. we only need to remove taxa from each random
    # tree. we do not need worry about taxa order as data.histo was already ordered in a way to match the order of the tree after its tips
    # are rid of unwanted taxa:
    drop.taxa <- tree.temp$tip.label[ ! tree.temp$tip.label %in% rownames(data.histo)]
    tree.temp <- drop.tip(tree.temp , drop.taxa)
    pfda <- phylo.fda.pred(XA,variable,tree.temp$tip.label,tree.temp,test.taxa,val = lambda, eqprior=TRUE, priin = 1)
    pfdas [[i]]<-pfda
    
    # ternary classification has two border lines. store the coefficients of the standard line equations:
    boundaryline1_coeffs[[i,"coeff_intercept"]]  <- pfda$fit$coefficients[[1,1]]
    boundaryline1_coeffs[[i,"coeff_compactness"]]  <- pfda$fit$coefficients[[2,1]]
    boundaryline1_coeffs[[i,"coeff_diameter"]]  <- pfda$fit$coefficients[[3,1]]
    
    boundaryline2_coeffs[[i,"coeff_intercept"]]  <- pfda$fit$coefficients[[1,2]]
    boundaryline2_coeffs[[i,"coeff_compactness"]]  <- pfda$fit$coefficients[[2,2]]
    boundaryline2_coeffs[[i,"coeff_diameter"]]  <- pfda$fit$coefficients[[3,2]]
    
    
    # some measures of performance for training from the confusion matrix
    rate.accuracy.terrestrial[i] <- pfda$confusion[1,1] / sum(pfda$confusion[,1])  # (truely predicted terrestrial)/(all terrestrial)
    rate.accuracy.amphibious[i]  <- pfda$confusion[2,2] / sum(pfda$confusion[,2])  # (truely predicted amphibious)/(all amphibious)
    rate.accuracy.aquatic[i]     <- pfda$confusion[3,3] / sum(pfda$confusion[,3])  # (truely predicted aquatic)/(all aquatic)
    rate.accuracy.overall[i]     <- sum(diag(pfda$confusion)) / sum(pfda$confusion)# (truely predicted instances) / (all instances),i.e. accuracy
    
    # 
    test.class <- as.matrix(predict(pfda, pfda$DATAtest, type="class"))
    predicted.ecologies[,i]<- test.class[,1]             # what are the predicted classes (0,1 or 2) of our test taxa using current tree
    
    # Discriminant scores scores for the fossils are retrieved this way.
    test.variates <- round(predict(pfda, pfda$DATAtest, type="variates"), digits = 2)
    
    # Posterior probabilities of group affiliations are available, too.
    test.prob <- round(predict(pfda, pfda$DATAtest, type="posterior"), digits = 2)
    probability.terrestrial[, i] <- test.prob[,1]
    probability.amphibious[, i] <- test.prob[,2]
    probability.aquatic[, i] <- test.prob[,3]
    
    # And now let's put everything together (here shown for the fossils):
    test.results <- cbind(test.class, test.prob, test.variates)
    colnames(test.results) <- c("predicted class", "P(terrestrial)", "P(amphibious)", "P(aquatic)", "DA1", "DA2")
    rownames(test.results) <- names(test.taxa)
    
    test.lambda.results[[i]]<- test.results
    
    # confusion matrix
    confusion.matrix[[i,1]] <- pfda$confusion[1,1]
    confusion.matrix[[i,2]] <- pfda$confusion[1,2]
    confusion.matrix[[i,3]] <- pfda$confusion[1,3]
    confusion.matrix[[i,4]] <- pfda$confusion[2,1]
    confusion.matrix[[i,5]] <- pfda$confusion[2,2]
    confusion.matrix[[i,6]] <- pfda$confusion[2,3]
    confusion.matrix[[i,7]] <- pfda$confusion[3,1]
    confusion.matrix[[i,8]] <- pfda$confusion[3,2]
    confusion.matrix[[i,9]] <- pfda$confusion[3,3]
  }      
  
  write.csv(boundaryline1_coeffs, file = "./BOUNDARY_LINE_1_COEFFS.csv")
  write.csv(boundaryline2_coeffs, file = "./BOUNDARY_LINE_2_COEFFS.csv")
  
  
  #number of times predicted terrestrial/amphibious/aquatic
  predictions.terrestrial <-row_count(as.data.frame(predicted.ecologies), count = "0")
  predictions.amphibious <-row_count(as.data.frame(predicted.ecologies), count = "1")
  predictions.aquatic <-row_count(as.data.frame(predicted.ecologies), count = "2")
  count.predictions.terrestrial <- predictions.terrestrial$rowcount ; names(count.predictions.terrestrial) <- names(test.taxa)
  count.predictions.amphibious <- predictions.amphibious$rowcount ; names(count.predictions.amphibious) <- names(test.taxa)
  count.predictions.aquatic <- predictions.aquatic$rowcount ; names(count.predictions.aquatic) <- names(test.taxa)
  
  # training performance metrics: median of accuracy rates for predicting each class label
  rates.accuracy <- cbind(median(rate.accuracy.overall), median(rate.accuracy.terrestrial), median(rate.accuracy.amphibious), median(rate.accuracy.aquatic)) 
  colnames(rates.accuracy) <- c("Training median overall accuracy", "Training median true terrestrial rate","Training median true amphibious rate", "Training median true aquatic rate")
  write.csv(rates.accuracy, file = "./TRAINING_PERFORMANCE.csv")
  
  #summarise median and maybe range values for predictions of all categories
  medians.terrestrial <- rowMedians(probability.terrestrial)
  medians.amphibious <- rowMedians(probability.amphibious)
  medians.aquatic <- rowMedians(probability.aquatic)
  medians.terrestrial          # median probability of terrestrial (over the given number of trees)
  medians.amphibious           # median probability of amphibious (over the given number of trees)
  medians.aquatic              # median probability of aquatic (over the given number of trees)
  
  
  write.csv(probability.terrestrial, file = "./TERRESTRIAL_PROBABILITIES.csv")
  write.csv(probability.amphibious, file = "./AMPHIBIOUS_PROBABILITIES.csv")
  write.csv(probability.aquatic, file = "./AQUATIC_PROBABILITIES.csv")
  
  
  summary.results <- matrix(ncol=8, nrow=length(test.set))
  for( i in 1:length(test.set) ) {
    all.median.probabilities <- c(medians.terrestrial[[test.set[[i]]]], medians.amphibious[[test.set[[i]]]], medians.aquatic[[test.set[[i]]]])
    predicted.class.index <- which.max( all.median.probabilities )     # which returns the index of max
    predicted.class <- colnames(pfda$confusion)[[predicted.class.index]]
    summary.results[[i, 1]] <- variable[[test.set[[i]]]]    # test taxa original labels 
    summary.results[[i, 2]] <- predicted.class                     # test taxa predicted labels
    summary.results[[i, 3]] <- count.predictions.terrestrial[[test.set[[i]]]]  # number of times test taxa are predicted to be terrestrial
    summary.results[[i, 4]] <- count.predictions.amphibious[[test.set[[i]]]]  # number of times test taxa are predicted to be amphibious
    summary.results[[i, 5]] <- count.predictions.aquatic[[test.set[[i]]]]  # number of times test taxa are predicted to be aquatic
    summary.results[[i, 6]] <- medians.terrestrial[[test.set[[i]]]] # median probability of terrestrial prediction
    summary.results[[i, 7]] <- medians.amphibious[[test.set[[i]]]] # median probability of amphibious prediction
    summary.results[[i, 8]] <- medians.aquatic[[test.set[[i]]]] # median probability of aquatic prediction
  }
  
  colnames(summary.results) <- c("actual class", "predicted class", "N(terrestrial)", "N(amphibious)", "N(aquatic)", "median P(terrestrial)", "median P(amphibious)", "median P(aquatic)")
  rownames(summary.results) <- test.set
  write.csv(summary.results, file = "./SUMMARY.csv")
  
  colnames(confusion.matrix) <- c("[1,1]", "[1,2]", "[1,3]", "[2,1]", "[2,2]", "[2,3]", "[3,1]", "[3,2]", "[3,3]")
  write.csv(confusion.matrix, file = "./CONFMAT.csv")
  
}


if (runStandAlone){
  result <- system.time(run_single())
  write.table(as.matrix(result), file=file.path(baseDir, "runtime.log"), row.names=T, col.names=F)
}






