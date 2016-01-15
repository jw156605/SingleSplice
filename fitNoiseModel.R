library(statmod)
suppressMessages(library(matrixStats))

fitMeanVarModel <- function(ERCC_rpkms_size_norm,cutoff)
{
	means = rowMeans(ERCC_rpkms_size_norm)
	vars = rowVars(ERCC_rpkms_size_norm)
	cv2 = vars/means^2
	useForFit = which(means > cutoff)
	fit = glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ), cv2[useForFit] )
	postscript("noise_model.eps")
	plot(means[which(means>0)],cv2[which(means>0)],log="xy",main="Mean-Variance Relationship",xlab="Mean Expression",ylab="CV^2")
	xg <- 10^seq( -1, 6, length.out=100 )
	lines( xg, coefficients(fit)["a0"] + coefficients(fit)["a1tilde"]/xg, col="red" )
	invisible(dev.off())
	return (c(coefficients(fit)["a0"],coefficients(fit)["a1tilde"]))
}

fitDropoutModel <- function(ERCC_rpkms_size_norm,cutoff)
{
	means = rowMeans(ERCC_rpkms_size_norm)
	num_dropouts = rowSums(ERCC_rpkms_size_norm==0)
	n = ncol(ERCC_rpkms_size_norm)
	total = rep(n,length(num_dropouts))
	df = data.frame(num_dropouts,total)
	fit = glm(cbind(num_dropouts,total-num_dropouts) ~ log(means+1),family=binomial(logit),data=df)
	b0 = fit$coefficients[1]
	b1 = fit$coefficients[2]
	p = num_dropouts/n
	sorted_means = sort(means)
	p_pred = (1/(1+exp(-(b0+b1*log(sorted_means+1)))))
	postscript("dropout_model.eps")
	plot(means+1,p,log="x",main="Dropout Rate",xlab="Mean Expression",ylab="Dropout Probability")
	lines(sorted_means+1,p_pred,col="red")
	invisible(dev.off())
	return (c(b0,b1))
}

args <- commandArgs(trailingOnly = TRUE)
ERCC_rpkms_size_norm = as.matrix(read.csv(args[1],header=TRUE,row.names=1))
cutoff = as.numeric(args[2])

meanVarModel = fitMeanVarModel(ERCC_rpkms_size_norm,cutoff)
dropoutModel = fitDropoutModel(ERCC_rpkms_size_norm)
cat("",meanVarModel[1],meanVarModel[2],dropoutModel[1],dropoutModel[2],"\n")