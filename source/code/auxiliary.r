# auxiliary.r
# These are some helper functions
# Cem Ozen, Aug. 3, 2022
library(stringr)

concord_tree_data <- function(tree, data){
  check <- name.check(tree, data.names=rownames(data))
  if (typeof(check) == "list"){ 
    dropList <- check$tree_not_data  
    if (length(dropList) > 0){
      tree <- drop.tip(tree, dropList)  # drop tips that do not have a match in data
      print("The following taxa are in the tree but not in the data. Hence we remove them from the tree:")
      print(dropList)
    }
    dropList <- check$data_not_tree
    if (length(dropList) > 0){
      data <- data[!(rownames(data) %in% dropList), ]
      print("The following taxa are in the data but not in the tree. Hence we remove them from the data:")
      print(dropList)
    }
  }  
  list("tree" = tree, "data" = data)
}


manipulateData  <- function(data, applyDataFilter){
  # applyDataFilter is any of: "Fabbri", "F0", F0D02 
  if (applyDataFilter == "Fabbri"){                   
    both <- data[, "finest.category"] == "both"
    data <- data[! both ,]
  } 
  else if (applyDataFilter == "Basic"){                   
    data[,"diving.or.not"] <- data[,"Diving"]
  }
  else if (applyDataFilter == "F0"){
    # we use Diving column in regression:
    data[,"diving.or.not"] <- data[,"Diving"]
    "keep only non-flying taxa"
    data <- data[data$Flying == 0,]
  }
  else if (applyDataFilter == "F0D02"){
    # we use Diving column in regression:
    data[,"diving.or.not"] <- data[,"Diving"]
    "keep only non-flying taxa"
    data <- data[data$Flying == 0,]
    "keep only D=0 and D=2 taxa. (retain UNKNOWN for test taxa)"
    data <- data[data$diving.or.not == 0 |  data$diving.or.not == 2 | data$diving.or.not == "UNKNOWN",]
  }
  else {
    stop("ERROR: manipulateData: undefined applyDataFilter is given!")
  }
  data
}


phyloDuplicateTip <- function(tree, species.label, new.species.label){
  # This function takes a tree and duplicates one of its tips. It labels
  # the new tip with the given name.
  # ----------------------------------------------------------------
  # Note that trees tips and nodes are labelled as such:
  # tips of a tree are enumerated from 1 to N_tip. Then the root of the tree is enumerated
  # N_tip+1 and the nodes are thus enumerated N_tip+1, N_tip+2, .., N_tip+N_node.
  # Another thing of importance is that for any ith edge, that is tree$edge[i,], there are
  # two components: tree$edge[i,1] gives the index of the node that is closer to the root of 
  # the tree and tree$edge[i,2] gives the index of the node that is further away from the 
  # root of the tree. If  tree$edge[i,2] <= tree$Nnode, due the enumeration convention
  # mentioned above the edge must have terminated at one of the tips of the tree. Otherwise
  # the edge is terminated in another node, so we must be still below (or left of the ) the branches of the tree.
  # if we view the tree from roots upwards (or from right to left).
  #-----------------------------------------------------------------
  # find index of a given tip by name. For example, find index of "t8"
  my.tip.index <- which(tree$tip.label == species.label)
  # find the node to which the tip is connected: node which "t8" is connected:
  # Note: tree$edge[i,2] stands for the ith edge's end node (among which a tip should be)
  #       and tree$edge[i,1] is the beginning node of the ith edge. 
  my.edge.index <- which(tree$edge[,2] == my.tip.index)
  my.edge.length <- tree$edge.length[my.edge.index]
  my.node <- tree$edge[my.edge.index,1]
  bind.tip(tree, new.species.label, edge.length = my.edge.length, where = my.node)
}


findSpeciesWithCommonRoot <- function(tree, species.label){
  # This function returns the names of other tips of a tree that are connected to the same
  # immediate node as the tip whose name is given. Note that it is possible that the given tip is connected
  # to a node that has its other forward edges connected to not tips but internal nodes. If that is the case,
  # it returns NA for internal nodes do not have a name. 
  my.tip.index <- which(tree$tip.label==species.label)
  # Find the index of the edge this tip is connected to:
  my.edge.index <- which(tree$edge[,2] == my.tip.index)
  # Find the node this edge is connected to:
  my.node <- tree$edge[my.edge.index,1]
  # Find all species (tips) connected to the same node:
  my.edges <- which(tree$edge[,1]==my.node)
  my.other.nodes <- tree$edge[my.edges,2]
  #my.other.nodes
  tree$tip.label[my.other.nodes]
}


get_base_duplicate_pairs <- function(data){
  # find all taxa that are duplicates of our normal taxa (assuming that they are "something.1" as in Spinosaurus_.1, Baryonyx.1 etc.)
  duplicates <- rownames(data[grep("\\.[0-9]+$", rownames(data)),])
  # extract their base names before the suffixes identifying them as duplicates:
  bases <- str_remove(duplicates, "\\.?[0-9]+$")
  # base_copy_pairs
  cbind(bases,duplicates)
}
