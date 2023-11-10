library(data.table)    


add_test_taxon_variant <- function(data, baseTaxon, suffix, MD, Log10MD, Cg, reference, notes){
  # WARNING: Never enter a comma in any string argument, as it would mess up the CVS file column structure
  variant <- data[rownames(data)==baseTaxon,]
  rownames(variant) <- paste0(baseTaxon,suffix)
  variant[,"MD..mm."] <- MD
  if ("MD.log" %in%  colnames(data)){   # they named ribs and femur data columns differently (darn it) so let's deal with it
    variant[,"MD.log"]  <- Log10MD    # for femur this work
  }else{
    variant[,"MD.Log"]  <- Log10MD    # for ribs this work
  } 
  variant[,"Global.compactness"]  <- Cg
  variant[,"Reference"] <- reference
  variant[,"Specimen"] <- "NA"
  variant[,"Type.of.data"] <- notes
  rbind(data, variant)
}

dataFile <- "c:/dino_nature/data/femur_compactness_all.csv"
newDataFile <-"c:/dino_nature/data/femur_DS1.csv"
data <- read.csv(dataFile, row.names = 1) 
dim(data)
data <- add_test_taxon_variant(data, "Spinosaurus_", ".1",  81.52,	1.91126,	0.804, "Nathan et al", "D same as Fabbri Cg low of 3 CT scans")
data <- add_test_taxon_variant(data, "Spinosaurus_", ".2",  81.52,	1.91126,	0.849, "Nathan et al", "D same as Fabbri Cg medium of 3 CT scans")
data <- add_test_taxon_variant(data, "Spinosaurus_", ".3",  81.52,	1.91126,	0.888, "Nathan et al", "D same as Fabbri Cg high of 3 CT scans")
data <- add_test_taxon_variant(data, "Spinosaurus_", ".4",  81.52,	1.91126,	0.914, "Nathan et al", "D same as Fabbri Cg Fabbri et al claim in reply")
data <- add_test_taxon_variant(data, "Spinosaurus_", ".5",  133.434,	2.12527,	0.968, "Nathan et al", "D scaled to adult Cg same as Fabbri")
data <- add_test_taxon_variant(data, "Spinosaurus_", ".6",  133.434,	2.12527,	0.804, "Nathan et al", "D scaled to adult Cg low of 3 CT scans")
data <- add_test_taxon_variant(data, "Spinosaurus_", ".7",  133.434,	2.12527,	0.849, "Nathan et al", "D scaled to adult Cg medium of 3 CT scans")
data <- add_test_taxon_variant(data, "Spinosaurus_", ".8",  133.434,	2.12527,	0.888, "Nathan et al", "D scaled to adult Cg high of 3 CT scans")
data <- add_test_taxon_variant(data, "Spinosaurus_", ".9",  133.434,	2.12527,	0.914, "Nathan et al", "D scaled to adult Cg Fabbri et al claim in reply")
data <- add_test_taxon_variant(data, "Spinosaurus_", ".10", 57.39,	1.75884, 0.699,    "Nathan et al", "D baby femur Cg SB scan high")
data <- add_test_taxon_variant(data, "Spinosaurus_", ".11", 57.39,	1.75884, 0.689,    "Nathan et al", "D baby femur Cg SB scan low")

data <- add_test_taxon_variant(data, "Suchomimus", ".1", 120.6,	2.08135,	0.628, "Nathan et al", "D same as Fabbri Cg from SB scan")
data <- add_test_taxon_variant(data, "Suchomimus", ".2", 146.4,	2.16554,	0.682, "Nathan et al", "D actual max GAD500 Cg same as Fabbri")
data <- add_test_taxon_variant(data, "Suchomimus", ".3", 146.4,	2.16554,	0.628, "Nathan et al", "D actual max GAD500 Cg from SB scan")


data <- add_test_taxon_variant(data, "Baryonyx", ".1", 154,	2.18752,	0.887,  "Nathan et al", "D same as Fabbri Cg scanned Cg high")
data <- add_test_taxon_variant(data, "Baryonyx", ".2", 154,	2.18752,	0.826, "Nathan et al", "D same as Fabbri Cg scanned Cg median")
data <- add_test_taxon_variant(data, "Baryonyx", ".3", 154,	2.18752,	0.767, "Nathan et al", "D same as Fabbri Cg scanned Cg low")

data[grep("Spinosaurus", rownames(data)),]
data[grep("Suchomimus", rownames(data)),]
data[grep("Baryonyx", rownames(data)),]

write.csv(data, file = newDataFile)

#************************************************************************************************************

dataFile <- "c:/dino_nature/data/ribs_compactness_all.csv"
newDataFile <-"c:/dino_nature/data/ribs_DS1.csv"
data <- read.csv(dataFile, row.names = 1) 
dim(data)
data <- add_test_taxon_variant(data, "Spinosaurus", ".1",  35.1,	1.54531,	0.8379,  "Nathan et al", "D same as Fabbri Cg scaled by 0.9")
data <- add_test_taxon_variant(data, "Spinosaurus", ".2",  57.4517,	1.7593,	0.931,   "Nathan et al", "D scaled to adult Cg same as Fabbri")
data <- add_test_taxon_variant(data, "Spinosaurus", ".3",  57.4517,	1.7593,	0.8379,  "Nathan et al", "D scaled to adult Cg scaled by 0.9")

data <- add_test_taxon_variant(data, "Baryonyx", ".1",  42.2,	1.62531,	0.8289,   "Nathan et al", "D same as Fabbri Cg scaled by 0.9")
write.csv(data, file = newDataFile)
