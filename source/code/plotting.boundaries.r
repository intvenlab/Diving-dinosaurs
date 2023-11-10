# plotting.boundaries.r
# after regression, we can plot decision boundary lines
# Cem Ozen, Aug 3, 2022.

# STILL WORKING ON THIS!
#*****************************************************

tree.no <- 1
test.compactness <- compactness[test.taxa]
test.diameter <- diameter[test.taxa]

# Boundary line coeffs:
coeff_intercept <- pfdas[[tree.no]]$fit$coefficients[[1]]
coeff_compactness <- pfdas[[tree.no]]$fit$coefficients[[2]]
coeff_diameter <- pfdas[[tree.no]]$fit$coefficients[[3]]

straightLineFunc <- function(x, a, b, c){y = -1/a * (b*x+c)}
decisionBoundaryFunc <- function(x) { straightLineFunc(x, coeff_compactness, coeff_diameter, coeff_intercept) }
decisionBoundaryLine <- lapply(training.diameter, decisionBoundaryFunc)

plot(x = training.diameter, y=training.compactness, pch=19, col=as.factor(diving.labels.T), xlab="diameter",ylab="compactness")
points(x=test.diameter, y=test.compactness, pch=2, col="blue", cex=1.0)
text(x=test.diameter+0.1, y=test.compactness, labels=names(test.taxa),cex=0.8)
lines(training.diameter, decisionBoundaryLine , type='l', col = "black", lwd=2,lty=1)


pfdas[[tree.no]]$confusion

# CEM: If we predict test taxa labels that we actually had to start with, the following can be used to evaluate the permormance of our prediction:
#confusion(as.factor(summary.results[,2]), as.factor(summary.results[,1]))  # confusion matrix
#summary.accuracy = sum(summary.results[,1] == summary.results[,2])/dim(summary.results)[[1]]
#summary.true.diving.rate = sum(summary.results[,2] == "2" & summary.results[,1] == "2") / sum(summary.results[,1] == "2")   # TP/(TP+FN) = (Pred = 2 AND Actual =2)/(Actual = 2)
#summary.true.nondiving.rate = sum(summary.results[,2] == "0" & summary.results[,1] == "0") / sum(summary.results[,1] == "0") # TN/(TN+FP) = (Pred = 0 AND Actual =0)/(Actual = 0)
#summary.accuracy
#summary.true.diving.rate
#summary.true.nondiving.rate
