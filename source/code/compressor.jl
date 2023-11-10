
include("misc.jl")
using .RunOrAnalyze

if (length(ARGS) != 2)
    error("Need 2 command-line arguments!")
end

runType = ARGS[1]             # any of: single, loocv, kfoldcv, bs
baseDir = ARGS[2]
setKey  = r"set_[0-9]+"       # relevant only for runs involving multiple sets (loocv, bootstrap, kfoldcv)

compressor(runType, baseDir)
