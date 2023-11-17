# Diving-dinosaurs
Diving dinosaurs? Caveats on the use of bone compactness and pFDA for inferring lifestyle.

Contents under revsion.

Contains a slightly modified version of phylo.fda.R:   https://github.com/lschmitz/phylo.fda

Install notes: 

Depends on:  python, r, r studio, julia language, jupyter notebook

Add support for the jupyter notebook julia kernel by installing it from Julia: 

```import Pkg; Pkg.add("IJulia"); Pkg.add("DataFrames"); Pkg.add("DataFrames"); Pkg.add("CSV"); Pkg.add("RCall"); Pkg.add("Glob"); ```

From R studio, install these additional packages via the IDE, under "Tools" -> "Install Packages"

```geomorph phytools geomorph geiger NBClust svgViewR formattable dispRity kader paleotree HDInterval mda qpcR klaR dplyr data.table lattice sjmisc lme4 modelr```

Path is assumed to be c:\dino_nature\code; if it's in a different location all source files must have the path replaced with the correct location. 

To execute: 

1.  Launch Jupyter notebook as follows: 

```
cd \dino_nature\code
jupyter notebook
```

2.  Pick the run_or_analyze.ipynb file

3.  Customize the job as desired.  It's a good idea to change your output directory after each run, as it will overwrite files if needed if run a second time on the same output directory resulting in unexpected structures.  

```

homeDir = "C:/dino_nature/output/2023-11-10-D4-Run2"    # where full analysis folders are (or will be) placed at

boneType            = "ribs"                                                                # enter "femur" or "ribs" 
dataFile            = "C:/dino_nature/data/D4_ribs.csv"                         # bone data file (with path). See above.
applyDataFilter     = "FALSE"                                                                # enter FALSE (no filter) or filter name. (see auxiliary.r) 
addTreeTips         = "TRUE"                                                                 # see explanations
Ntrees              = 100                                                                    # no of random trees to be generated
calculateTrees      = "TRUE"                                                                 # calculate (TRUE) or read (FALSE) trees
saveTrees           = "FALSE"                                                                # save (TRUE) trees (if they are calculated here)
Kcrossvalidation    = 10                                                                     # # of data partitions (K) in a K-fold CV 
Nkfoldcv            = 10                                                                     # how many rounds of k-fold cv is wanted   
Nbootstrapsets      = 2000                                                                   # # of bootstrap sets
```






