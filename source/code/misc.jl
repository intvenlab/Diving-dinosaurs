
module RunOrAnalyze
#=
this module provides the function `compressor` which compresses PROBABILITIES.csv and CONFMAT.csv files created by single prediction, bootstrap, loocv, and kfoldcv type runs 
involving binary classification. If ternary classification is done, the module has to 
be generalzied for that task.

Also it performs full analysis, which  means performing single (binary_classification), loocv, k-fold cv, and bootstrap runs 
altogether called full suit.

Explanation for when addTreeTips is to be set TRUE. if
If either holds: 
                              a) user added a test taxon as a variant of an existing taxon 
                              b) any bootstrap run. 
in a) user may want to change a value of one of the predictor fields of a certain taxon to see how such a change affect
the prediction of the diving label. (For example, data may have Spinosaurus_ and we may want to test also a version of Spinosaurus_ which
we can add to data by changing compactness or diameter fields of Spinosaurus_. Such a line need to be added as in Spinosaurus_.1 This new line
in data needs to be addressed to Spinosaurus_ in the phylo tree. We do this by artificially creating a tip named Spinosaurus_.1 in the tree that is 
connected to the same root as Spinosaurus_ with same edge length.)
Similarly, in b) same necessity emerges. This time whether such a user added taxon exists or not, the training part of the data will be resampled with replacement.
meaning that there will be repeated taxa in the data. These taxons are named by R in the fashion described above by adding a dot and an integer to the 
taxons name with which it is identical. These type of copycat taxa need their own tips in the tree which are to be copycat tips of the copied taxon's tip.
Therefore, each time we do bootstrap, we automatically addTreeTips to TRUE (see .RunOrAnalayze.bootstrap function)   

Always use the designated data directory for data file location. Easier to keep track of stuff. 

# Explanation for the functions: 
#     single:
#     loocv:
#     kfoldcv:
#     bootstrap:

Cem Ozen, Aug. 11, 2022.

=#
using DataFrames, CSV, RCall
export compressor, single, loocv, kfoldcv, bootstrap, newfullsuitrun, fullsuitpostprocess

# Functions for running or analyzing. 

function single(homeDir, boneType, dataFile, applyDataFilter, addTreeTips, Ntrees, calculateTrees,saveTrees)
    # start an R session in the background, load the R code containing code and input 
    R"source(\"c:/dino_nature/code/binary_classifier.r\")"
    # check subdirectory "single", if not exists, create 
    baseDir = joinpath(homeDir, "single")
    isdir(baseDir) || mkdir(baseDir)
    # update input to override those in binary_classifier.r
    @rput baseDir boneType dataFile applyDataFilter addTreeTips Ntrees calculateTrees saveTrees
    # run the function 
    R"run_single()"
end

function loocv(homeDir, boneType, dataFile, applyDataFilter, addTreeTips, Ntrees, calculateTrees,saveTrees)
    # start an R session in the background, load the R code containing code and input 
    R"source(\"c:/dino_nature/code/loocv_binary.r\")"   
    # check subdirectory "loocv", if not exists, create
    baseDir = joinpath(homeDir, "loocv")
    isdir(baseDir) || mkdir(baseDir)
    # update input to override those in loocv_binary.r
    @rput baseDir boneType dataFile applyDataFilter addTreeTips Ntrees calculateTrees saveTrees
    # run the function
    R"run_loocv()"
end

function kfoldcv(homeDir, boneType, dataFile, applyDataFilter, addTreeTips, Ntrees, calculateTrees,saveTrees, Kcrossvalidation, Nkfoldcv, ikfoldcv)
    # start an R session in the background, load the R code containing code and input 
    R"source(\"c:/dino_nature/code/kfoldcv_binary.r\")"  
    # check subdirectories "kfoldcv" or "kfoldcv_1" etc, if not exists, create
    if Nkfoldcv == 1
        baseDir = joinpath(homeDir, "kfoldcv") 
        isdir(baseDir) || mkdir(baseDir)
        # update input to override those in kfoldcv_binary.r
        @rput baseDir boneType dataFile applyDataFilter addTreeTips Ntrees calculateTrees saveTrees Kcrossvalidation
        # run the function
        R"run_kfoldcv()"
    else         # Nkfoldcv > 1 is implied
        baseDir = joinpath(homeDir, "kfoldcv_$ikfoldcv") 
        isdir(baseDir) || mkdir(baseDir)
        # update input to override those in kfoldcv_binary.r
        @rput baseDir boneType dataFile applyDataFilter addTreeTips Ntrees calculateTrees saveTrees Kcrossvalidation
        # run the function
        R"run_kfoldcv()"
    end 
end

function bootstrap_sets_generator(baseDir, boneType, dataFile, applyDataFilter, Nbootstrapsets)
    # start an R session in the background, load the R code containing code and input
    R"source(\"c:/dino_nature/code/bootstrap_sets_generator.r\")"   
    # update input to override those in bootstrap_sets_generator.r
    @rput baseDir boneType dataFile applyDataFilter Nbootstrapsets
    # run the function
    R"create_bootstrap_sets()"
end

function bootstrap(homeDir, boneType, dataFile, applyDataFilter, Ntrees, calculateTrees,saveTrees, Nbootstrapsets)
    addTreeTips = "TRUE"     # always the case for bootstrap (see explanations above)
    # start an R session in the background, load the R code containing code and input
    R"source(\"c:/dino_nature/code/bootstrap_binary.r\")" 
    # check subdirectory "bs", if not exists, create
    baseDir = joinpath(homeDir, "bs") 
    isdir(baseDir) || mkdir(baseDir)
    # create bootstrap sets within "bs" subdirectory
    bootstrap_sets_generator(baseDir, boneType, dataFile, applyDataFilter, Nbootstrapsets)
    # update input to override those in bootstrap_binary.r
    @rput baseDir boneType applyDataFilter addTreeTips Ntrees calculateTrees saveTrees
    # run the function
    R"run_bootstrap()"
end


function newfullsuitrun(homeDir, boneType, dataFile, applyDataFilter, addTreeTips, Ntrees, calculateTrees,saveTrees,Kcrossvalidation,Nkfoldcv,Nbootstrapsets)
    # check if directory to host the full statistical analysis exists, if not create
    isdir(homeDir) || mkdir(homeDir)
    # perform full statistical analysis: single, loocv, k-foldcv, bootstrap 
    single(homeDir, boneType, dataFile, applyDataFilter, addTreeTips, Ntrees, calculateTrees,saveTrees)
    loocv(homeDir, boneType, dataFile, applyDataFilter, addTreeTips, Ntrees, calculateTrees,saveTrees)
    for ikfoldcv ∈ 1:Nkfoldcv
        kfoldcv(homeDir, boneType, dataFile, applyDataFilter, addTreeTips, Ntrees, calculateTrees,saveTrees, Kcrossvalidation,Nkfoldcv, ikfoldcv)
    end
    bootstrap(homeDir, boneType, dataFile, applyDataFilter, Ntrees, calculateTrees,saveTrees, Nbootstrapsets)
end

function fullsuitpostprocess(homeDir, Nkfoldcv)
    isdir(homeDir) || return "directory does not exist!"
    # create post process results folders
    outputDir = joinpath(homeDir, "analysis") 
    isdir(outputDir) || mkdir(outputDir)
    # single prediction 
    baseDir = joinpath(homeDir, "single")
    isdir(baseDir) && compressor("single", baseDir)
    isdir(baseDir) && cp(joinpath(baseDir,"singlepredictionprobabilites.csv"), joinpath(outputDir,"singlepredictionprobabilites.csv"), force=true)
    isdir(baseDir) && cp(joinpath(baseDir,"singlepredictionconfusionmatrices.csv"), joinpath(outputDir,"singlepredictionconfusionmatrices.csv"), force=true)
    isdir(baseDir) && cp(joinpath(baseDir,"singlePCM.csv"), joinpath(outputDir,"singlePCM.csv"), force=true)
    # loocv
    baseDir = joinpath(homeDir, "loocv")
    isdir(baseDir) && compressor("loocv", baseDir)
    isdir(baseDir) && cp(joinpath(baseDir,"loocvpredictionprobabilites.csv"), joinpath(outputDir,"loocvpredictionprobabilites.csv"), force=true)
    isdir(baseDir) && cp(joinpath(baseDir,"loocvpredictionconfusionmatrices.csv"), joinpath(outputDir,"loocvpredictionconfusionmatrices.csv"), force=true)
    isdir(baseDir) && cp(joinpath(baseDir,"loocvPCM.csv"), joinpath(outputDir,"loocvPCM.csv"), force=true)
    # kfoldcv (Nkfoldcv times)
    if Nkfoldcv == 1
        baseDir = joinpath(homeDir, "kfoldcv")
        isdir(baseDir) && compressor("kfoldcv", baseDir)
        isdir(baseDir) && cp(joinpath(baseDir,"kfoldcvpredictionprobabilites.csv"), joinpath(outputDir,"kfoldcvpredictionprobabilites.csv"), force=true)
        isdir(baseDir) && cp(joinpath(baseDir,"kfoldcvpredictionconfusionmatrices.csv"), joinpath(outputDir,"kfoldcvpredictionconfusionmatrices.csv"), force=true)
        isdir(baseDir) && cp(joinpath(baseDir,"kfoldcvPCM.csv"), joinpath(outputDir,"kfoldcvPCM.csv"), force=true)
    else
        for i ∈ 1:Nkfoldcv
            baseDir = joinpath(homeDir, "kfoldcv_$i")
            isdir(baseDir) && compressor("kfoldcv", baseDir)
            isdir(baseDir) && cp(joinpath(baseDir,"kfoldcvpredictionprobabilites.csv"), joinpath(outputDir,"kfoldcvpredictionprobabilites_$i.csv"), force=true)
            isdir(baseDir) && cp(joinpath(baseDir,"kfoldcvpredictionconfusionmatrices.csv"), joinpath(outputDir,"kfoldcvpredictionconfusionmatrices_$i.csv"), force=true)
            isdir(baseDir) && cp(joinpath(baseDir,"kfoldcvPCM.csv"), joinpath(outputDir,"kfoldcvPCM_$i.csv"), force=true)
        end
    end
    # bootstrap
    baseDir = joinpath(homeDir, "bs")
    isdir(baseDir) && compressor("bs", baseDir)
    isdir(baseDir) && cp(joinpath(baseDir,"bspredictionprobabilites.csv"), joinpath(outputDir,"bspredictionprobabilites.csv"), force=true)
    isdir(baseDir) && cp(joinpath(baseDir,"bspredictionconfusionmatrices.csv"), joinpath(outputDir,"bspredictionconfusionmatrices.csv"), force=true)
    isdir(baseDir) && cp(joinpath(baseDir,"bsPCM.csv"), joinpath(outputDir,"bsPCM.csv"), force=true)
end



# Functions for compressing output 

searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

function compressRow(row)
    tmp = [(count(==(x), row),x) for x in Set(row)]
end

function compress_probabilities(data)
    fileStr = ""
    for i ∈ 1:nrow(data)
        fileStr = fileStr * data[i,1] * ", " * string(compressRow(data[i,2:end])) * "\n"
    end
    fileStr = replace(fileStr, "[" => "{")
    fileStr = replace(fileStr, "]" => "}")
    fileStr = replace(fileStr, "(" => "{")
    fileStr = replace(fileStr, ")" => "}")
    return fileStr
end

function compress_confusionmatrix(data)
    V = Vector([Vector(data[i,2:5]) for i in 1:nrow(data)])
    fileStr = string(compressRow(V)) * "\n"
    fileStr = replace(fileStr, "[" => "{")
    fileStr = replace(fileStr, "]" => "}")
    fileStr = replace(fileStr, "(" => "{")
    fileStr = replace(fileStr, ")" => "}")
    return fileStr
end

function compress_confmat_and_prob(dataP, dataC)
    Ntaxa = nrow(dataP)     # number of test taxa, from the probability file
    Ntrees = ncol(dataP)-1    # number of trees from the probability file
    check1 = (nrow(dataC) == Ntrees)   # number of trees from the confusion matrix should match
    check2 = all(all.( x->typeof(x)<:Integer  ,eachcol(dataC) ))  # is the confusion matrix made of integers
    check3 = all(all.( x->typeof(x)<:Real  ,eachcol(dataP[:,2:end]) ))  # are the probabilities all Real?
    if check1 && check2 && check3
        fileStr = ""
        pcm = Array{Any}(undef, Ntaxa, Ntrees)
        for taxon ∈ 1:Ntaxa
            for tree ∈ 1:Ntrees
                pcm[taxon, tree] = (Array(dataC[tree,2:5]), dataP[taxon,tree+1])
            end
            fileStr = fileStr * dataP[taxon,1] * ", " * string(compressRow(pcm[taxon,:])) * "\n"
        end
        fileStr = replace(fileStr, "[" => "{")
        fileStr = replace(fileStr, "]" => "}")
        fileStr = replace(fileStr, "(" => "{")
        fileStr = replace(fileStr, ")" => "}")
        return fileStr
    else
        return ""     # if confusion matrix or probabilities contain non-numeric entries
    end
end

function compressor(runType, baseDir)
    if runType == "single"    # single prediction mode
        outFileName1 = "singlepredictionprobabilites.csv"
        outFileName2 = "singlepredictionconfusionmatrices.csv"
        outFileName3 = "singlePCM.csv"
    elseif runType == "bs"    # bootstrap mode
        outFileName1 = "bspredictionprobabilites.csv"
        outFileName2 = "bspredictionconfusionmatrices.csv"
        outFileName3 = "bsPCM.csv"
    elseif runType == "loocv"    # bootstrap mode
        outFileName1 = "loocvpredictionprobabilites.csv"
        outFileName2 = "loocvpredictionconfusionmatrices.csv"
        outFileName3 = "loocvPCM.csv"
    elseif runType == "kfoldcv"    # bootstrap mode
        outFileName1 = "kfoldcvpredictionprobabilites.csv"
        outFileName2 = "kfoldcvpredictionconfusionmatrices.csv"
        outFileName3 = "kfoldcvPCM.csv"
    else 
        error("runType is not defined!")
    end
    cd(baseDir)
    if runType == "single"
        # input and output files with their paths:
        inpFile = joinpath(baseDir, "PROBABILITIES.csv")
        dataP = DataFrame(CSV.File(inpFile))

        outFile = joinpath(baseDir, outFileName1)
        fileStr = compress_probabilities(dataP)
        write(outFile, fileStr)     # writing compressed probability output
    
        inpFile = joinpath(baseDir, "CONFMAT.csv")
        dataC = DataFrame(CSV.File(inpFile))

        outFile = joinpath(baseDir, outFileName2)
        fileStr = compress_confusionmatrix(dataC)
        write(outFile, fileStr)     # writing compressed confusion matrix output

        outFile = joinpath(baseDir, outFileName3)
        fileStr = compress_confmat_and_prob(dataP, dataC)
        write(outFile, fileStr)     # writing compressed combinations of probability and confusion matrix 


    else           # this includes bs, loocv, kfoldcv type runs involving sets
        setKey  = r"set_[0-9]+"       
        setList = searchdir(baseDir, setKey)
        outFile = joinpath(baseDir, outFileName1)
        fileStr = ""
        for dir ∈ setList
            inpFile = joinpath(baseDir, dir, "PROBABILITIES.csv")
            if !isfile(inpFile)
                println(inpFile * "does not exist! Skipping!\n")
                continue
            end
            data = DataFrame(CSV.File(inpFile))
            fileStr = fileStr * compress_probabilities(data) 
        end
        write(outFile, fileStr)

        outFile = joinpath(baseDir, outFileName2)
        fileStr = ""
        for dir ∈ setList
            inpFile = joinpath(baseDir, dir, "CONFMAT.csv")
            if !isfile(inpFile)
                println(inpFile * "does not exist! Skipping!\n")
                continue
            end
            data = DataFrame(CSV.File(inpFile))
            fileStr = fileStr * compress_confusionmatrix(data) 
        end
        write(outFile, fileStr)
        
        outFile = joinpath(baseDir, outFileName3)
        fileStr = ""
        for dir ∈ setList
            inpFileP = joinpath(baseDir, dir, "PROBABILITIES.csv")
            inpFileC = joinpath(baseDir, dir, "CONFMAT.csv")
            if !isfile(inpFileP) || !isfile(inpFileC)
                println("Either probability or confusion matrix is missing! Skipping a set!\n")
                continue
            end
            dataP = DataFrame(CSV.File(inpFileP))
            dataC = DataFrame(CSV.File(inpFileC))
            fileStr = fileStr * compress_confmat_and_prob(dataP, dataC)
        end
        write(outFile, fileStr)
    end
end

end












