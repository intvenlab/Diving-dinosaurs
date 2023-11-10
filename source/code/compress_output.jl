#=
compress_output.jl

this code works for single prediction, bootstrap, loocv, and kfoldcv type runs 
involving binary classification. If ternary classification is done, this codehas to 
be generalzied for that task.

Cem Ozen, Aug. 11, 2022.
=#
using DataFrames, CSV

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

function compressor(runType, baseDir)
    if runType == "single"    # single prediction mode
        outFileName1 = "singlepredictionprobabilites.csv"
        outFileName2 = "singlepredictionconfusionmatrices.csv"
    elseif runType == "bs"    # bootstrap mode
        outFileName1 = "bspredictionprobabilites.csv"
        outFileName2 = "bspredictionconfusionmatrices.csv"
    elseif runType == "loocv"    # bootstrap mode
        outFileName1 = "loocvpredictionprobabilites.csv"
        outFileName2 = "loocvpredictionconfusionmatrices.csv"
    elseif runType == "kfoldcv"    # bootstrap mode
        outFileName1 = "kfoldcvpredictionprobabilites.csv"
        outFileName2 = "kfoldcvpredictionconfusionmatrices.csv"
    else 
        error("runType is not defined!")
    end
    cd(baseDir)
    if runType == "single"
        # input and output files with their paths:
        inpFile = joinpath(baseDir, "PROBABILITIES.csv")
        outFile = joinpath(baseDir, outFileName1)
        data = DataFrame(CSV.File(inpFile))
        fileStr = compress_probabilities(data)
        write(outFile, fileStr)     # writing compressed probability file
    
        inpFile = joinpath(baseDir, "CONFMAT.csv")
        outFile = joinpath(baseDir, outFileName2)
        data = DataFrame(CSV.File(inpFile))
        fileStr = compress_confusionmatrix(data)
        write(outFile, fileStr)
    else           # this includes bs, loocv, kfoldcv type runs involving sets
        setList = searchdir(baseDir, setKey)
        outFile = joinpath(baseDir, outFileName1)
        fileStr = ""
        for dir ∈ setList
            inpFile = joinpath(baseDir, dir, "PROBABILITIES.csv")
            data = DataFrame(CSV.File(inpFile))
            fileStr = fileStr * compress_probabilities(data) 
        end
        write(outFile, fileStr)
        outFile = joinpath(baseDir, outFileName2)
        fileStr = ""
        for dir ∈ setList
            inpFile = joinpath(baseDir, dir, "CONFMAT.csv")
            data = DataFrame(CSV.File(inpFile))
            fileStr = fileStr * compress_confusionmatrix(data) 
        end
        write(outFile, fileStr)
    end
end


if (length(ARGS) != 2)
    error("Need 2 command-line arguments!")
end

runType = ARGS[1]             # any of: single, loocv, kfoldcv, bs
baseDir = ARGS[2]
setKey  = r"set_[0-9]+"       # relevant only for runs involving multiple sets (loocv, bootstrap, kfoldcv)

compressor(runType, baseDir)








