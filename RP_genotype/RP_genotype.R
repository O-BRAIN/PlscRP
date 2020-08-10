# ______________________________________________________________________
# Using PLSC-RP for optimized analysis of high-dimensional genotype data
# ______________________________________________________________________


# Clear the environment
rm(list=ls())
library(far)
library(proxy)


###################################################################################################################################
# insert filename of example data here


# Example Data: 
filenameX1 <- "/.../X1.RData"
filenameX2 <- "/.../X2.RData"
filename_X1names <- "/.../X1names.RData"
filename_X2names <- "/.../X2names.RData"


# insert filename of results here
filenameRes <- "/.../PlscRP/RP_genotype/ExampleData/"


# orthogonalization of random matrix using Gram-Schmidt algorithm
orthonorm <- FALSE


# number of permutations
nperms <- 5000
# more accurate results for high number of permutations
# faster results for small number of permutations


# set working directory
setwd("/.../PlscRP/functions")

###################################################################################################################################



###################################################################################################################################
# Start of R Script
###################################################################################################################################

#-----------------------------
# Read in necessary  functions
#-----------------------------

source("Menon_lowerBound.R")
source("PLSC.R")
source("Out_of_sample_cov.R")
source("Similarity.R")



#------------------
# Load example data
#------------------

load(filenameX1)               # body height and serum vaspin concentration (phenotype data)
load(filenameX2)               # SNP data (genotype data)
load(filename_X1names)         # height, vaspin, ln(vaspin)
load(filename_X2names)         # SNP names



#-------------------------------------
# Check dimensionality of example data
#-------------------------------------

N <- nrow(X1)
D1 <- ncol(X1)
D2 <- ncol(X2)

print("Dimensionality of example data")
print("==============================")
print(paste("Sample size:", N))
print(paste("Number of phenotype measures:", D1))
print(paste("Number of SNPs:", D2))



# --------------------------
# Columnwise standardization
# --------------------------

X1.zscore <- scale(X1, center=TRUE, scale=TRUE)
X2.zscore <- scale(X2, center=TRUE, scale=TRUE)



# --------------------
# Cross-product matrix
# --------------------

rX1X2 <- cor(X1.zscore, X2.zscore)



# -------------------------------------
# RP for dimensionality reduction in X2
# -------------------------------------

# define dimensionality of low-dimensional subspace using Menon theorem
epsi <- 1.0    # distance preservation factor
prob <- 0.05   # 1 - probability of distance preservation

kMen <- lbMenon(N, epsi, prob)
kMen <- kMen$k

if(ceiling(kMen) < N){
  K <- N
}else{
  K <- kMen
}

print("Menon theorem")
print("=============")
print(paste("Menon lower bound:", kMen))
print(paste("Specified dimensionality of low-dimensional subspace:", K))


# create i.i.d. N(0,1) matrix
randmat <- matrix(rnorm(D2*K, mean=0, sd=1), D2, K)   # rows and columns are N(0,1)-distributed

# orthogonalization using Gram-Schmidt algorith
if(orthonorm){
  randmat <- orthonormalization(randmat, basis=FALSE, norm=TRUE)
}

# dimension reduction using RP
X2.RP <- X2.zscore %*% randmat



#------------
# Run PLSC-RP
#------------

# Permutation test
comp <- 1        # component for which permutation test should be performed

perm.PLSC_RP <- PLSC.permute(X1.zscore, X2.RP, nperms, comp)

print("Permutation test for PLSC-RP")
print("============================")
print(paste("Covariance of latent variables:", perm.PLSC_RP$cov))
print(paste("p-value of covariance:", perm.PLSC_RP$pval_cov))


# PLSC-RP analysis
out.PLSC_RP <- PLSC(X1.zscore, X2.RP)

w_X1 <- out.PLSC_RP$A
w_X2.RP <- out.PLSC_RP$B
LVX1 <- out.PLSC_RP$LV_X
LVX2.RP <- out.PLSC_RP$LV_Y
covval <- out.PLSC_RP$d



#--------------------------------------------------------------
# Select as many components to explain at least 80% of variance
#--------------------------------------------------------------

totalVar <- out.PLSC_RP$totalVar

varsum <- numeric(length(totalVar))
for(i in 1:length(totalVar)){
  varsum[i] <- sum(totalVar[1:i]) 
}
sel <- which(varsum >= 80)
sel <- sel[1]

w_X1 <- w_X1[,1:sel]
w_X2.RP <- w_X2.RP[,1:sel]
LVX1 <- LVX1[,1:sel]
LVX2.RP <- LVX2.RP[,1:sel]
covval <- covval[1:sel]

print("Results of PLSC-RP")
print("==================")
print(paste("Number of components:", sel))
print(paste("Variance explained:", sum(totalVar[1:sel])))
print(paste("Covariance of latent variables:", covval))



#------------------------------------------
# Out-of-sample covariance using 10-fold CV
#------------------------------------------

cov_oos_PLSC_RP <- Cov_oos_CV(X1.zscore, X2.RP, sel)

cov_oos_mean <- cov_oos_PLSC_RP$cov_oos_mean
causalComp <- which.max(diff(cov_oos_mean))
raw_cov_oos <- diff(cov_oos_mean)

print("Out-of-sample covariance using 10-fold CV")
print("=========================================")
print(paste("Causal Component:", causalComp))
print(paste("Out-of-sample covariance of causal component:", raw_cov_oos[causalComp]))



#-----------------------------------------------
# Backtransformation to get original SNP weights
#-----------------------------------------------

w_X1 <- w_X1[,causalComp]
w_X1 <- matrix(w_X1)
w_X2 <- t(rX1X2) %*% w_X1 %*% (1/(sum(w_X1^2) * covval[causalComp]))
colnorm <- sqrt(sum(w_X2^2))
w_X2 <- (1/colnorm) * w_X2
LVX2 <- LVX2.RP[,causalComp]



#--------------------------
# Saving results of PLSC-RP
#--------------------------
save(w_X1, file = paste(filenameRes, "w_X1_PLSC_RP.RData", sep=""))
save(w_X2, file = paste(filenameRes, "w_X2_PLSC_RP.RData", sep=""))
save(LVX1, file = paste(filenameRes, "LVX1_PLSC_RP.RData", sep=""))
save(LVX2, file = paste(filenameRes, "LVX2_PLSC_RP.RData", sep=""))



#--------------------------
# Presenting results for X1
#--------------------------

t <- 1:nrow(w_X1)

# causal phenotypes
maxval <- max(abs(w_X1[,causalComp]))     # maximum weight
highval <- maxval * 0.8                   # top 80%
highpheno <- which(abs(w_X1[,causalComp]) > highval)

print("Causal phenotypes by PLSC-RP")
print("============================")
print(paste("causal phenotypes:", X1names[highpheno]))

graphics.off()
X11()

dev.copy(pdf, file = paste(filenameRes, "w_X1_PLSC_RP.pdf", sep=""))

par(mar=c(5.1,4.5,4.1,2.1))
plot(t, w_X1[,causalComp], type="h", col="darkblue", ann = FALSE, cex.axis=2)
title(xlab="phenotype", cex.lab = 2)
title(ylab="phenotype weight", cex.lab = 2)
lines(t[highpheno], w_X1[highpheno,causalComp], type="h", col="red", lwd=1.5)
points(t, w_X1[,causalComp], pch=16, col="darkblue")
points(t[highpheno], w_X1[highpheno,causalComp], pch=16, col="red")
legend("bottomleft", X1names[highpheno], col = "red", cex = 1, lty = 1, lwd=c(1.5))
title(main="Phenotype weight profile of PLSC-RP", cex.main = 2)

dev.off()



#--------------------------
# Presenting results for X2
#--------------------------

t <- 1:nrow(w_X2)

# causal SNPs
maxval <- max(abs(w_X2[,causalComp]))     # maximum weight
highval <- maxval * 0.8                   # top 80%
highsnp <- which(abs(w_X2[,causalComp]) > highval)

print("Causal SNPs by PLSC-RP")
print("======================")
for(i in 1:length(highsnp)){
  print(paste("ID of causal SNPs: ", highsnp[i], ", name of causal SNPs: ", X2names[i], sep=""))
}

# zoomed in SNP weight profile presented in the manuscript
graphics.off()
X11()

dev.copy(pdf, file = paste(filenameRes, "w_X2_PLSC_RP.pdf", sep=""))

par(mar=c(5.1,4.5,4.1,2.1))
plot(t, w_X2[,causalComp], type="h", col="darkblue", ann = FALSE, cex.axis=2)
title(xlab="SNP", cex.lab = 2)
title(ylab="SNP weight", cex.lab = 2)
lines(t[highsnp], w_X2[highsnp,causalComp], type="h", col="red", lwd=1.5)
points(t, w_X2[,causalComp], pch=16, col="darkblue")
points(t[highsnp], w_X2[highsnp,causalComp], pch=16, col="red")
legend("bottomright", X2names, col = c("red"), cex = 1, lty = 1, lwd=1.5)
title(main="SNP weight profile of PLSC-RP", cex.main = 2)

dev.off()



#------------------------------------
# Run traditional PLSC for comparison
#------------------------------------

# Permutation test
comp <- 1        # component for which permutation test should be performed

perm.PLSC <- PLSC.permute(X1.zscore, X2.zscore, nperms, comp)

print("Permutation test for PLSC")
print("=========================")
print(paste("Covariance of latent variables:", perm.PLSC$cov))
print(paste("p-value of covariance:", perm.PLSC$pval_cov))


# PLSC analysis
out.PLSC <- PLSC(X1.zscore, X2.zscore)

w_X1 <- out.PLSC$A
w_X2 <- out.PLSC$B
LVX1 <- out.PLSC$LV_X
LVX2 <- out.PLSC$LV_Y
covval <- out.PLSC$d



#--------------------------------------------------------------
# Select as many components to explain at least 80% of variance
#--------------------------------------------------------------

totalVar <- out.PLSC$totalVar

varsum <- numeric(length(totalVar))
for(i in 1:length(totalVar)){
  varsum[i] <- sum(totalVar[1:i]) 
}
sel <- which(varsum >= 80)
sel <- sel[1]

w_X1 <- w_X1[,1:sel]
w_X2 <- w_X2[,1:sel]
LVX1 <- LVX1[,1:sel]
LVX2 <- LVX2[,1:sel]
covval <- covval[1:sel]


print("Results of PLSC")
print("===============")
print(paste("Number of components:", sel))
print(paste("Variance explained:", sum(totalVar[1:sel])))
print(paste("Covariance of latent variables:", covval))



#------------------------------------------
# Out-of-sample covariance using 10-fold CV
#------------------------------------------

cov_oos_PLSC <- Cov_oos_CV(X1.zscore, X2.zscore, sel)

cov_oos_mean <- cov_oos_PLSC$cov_oos_mean
causalComp <- which.max(diff(cov_oos_mean))
raw_cov_oos <- diff(cov_oos_mean)

print("Out-of-sample covariance using 10-fold CV")
print("=========================================")
print(paste("Causal Component:", causalComp))
print(paste("Out-of-sample covariance of causal component:", raw_cov_oos[causalComp]))



#--------------------------
# Saving results of PLSC-RP
#--------------------------

save(w_X1, file = paste(filenameRes, "w_X1_PLSC.RData", sep=""))
save(w_X2, file = paste(filenameRes, "w_X2_PLSC.RData", sep=""))
save(LVX1, file = paste(filenameRes, "LVX1_PLSC.RData", sep=""))
save(LVX2, file = paste(filenameRes, "LVX2_PLSC.RData", sep=""))



#--------------------------
# Presenting results for X1
#--------------------------

t <- 1:nrow(w_X1)

# causal phenotypes
maxval <- max(abs(w_X1[,causalComp]))     # maximum weight
highval <- maxval * 0.8                   # top 80%
highpheno <- which(abs(w_X1[,causalComp]) > highval)

print("Causal phenotypes by PLSC")
print("=========================")
print(paste("causal phenotypes:", X1names[highpheno]))

graphics.off()
X11()

dev.copy(pdf, file = paste(filenameRes, "w_X1_PLSC.pdf", sep=""))

par(mar=c(5.1,4.5,4.1,2.1))
plot(t, w_X1[,causalComp], type="h", col="darkblue", ann = FALSE, cex.axis=2)
title(xlab="phenotype", cex.lab = 2)
title(ylab="phenotype weight", cex.lab = 2)
lines(t[highpheno], w_X1[highpheno,causalComp], type="h", col="red", lwd=1.5)
points(t, w_X1[,causalComp], pch=16, col="darkblue")
points(t[highpheno], w_X1[highpheno,causalComp], pch=16, col="red")
legend("bottomleft", X1names[highpheno], col = "red", cex = 1, lty = 1, lwd=c(1.5))
title(main="Phenotype weight profile of PLSC", cex.main = 2)

dev.off()



#--------------------------
# Presenting results for X2
#--------------------------

t <- 1:nrow(w_X2)

# causal SNPs
maxval <- max(abs(w_X2[,causalComp]))     # maximum weight
highval <- maxval * 0.8                   # top 80%
highsnp <- which(abs(w_X2[,causalComp]) > highval)

print("Causal SNPs by PLSC")
print("===================")
for(i in 1:length(highsnp)){
  print(paste("ID of causal SNPs: ", highsnp[i], ", name of causal SNPs: ", X2names[i], sep=""))
}

# zoomed in SNP weight profile presented in the manuscript
graphics.off()
X11()

dev.copy(pdf, file = paste(filenameRes, "w_X2_PLSC.pdf", sep=""))

par(mar=c(5.1,4.5,4.1,2.1))
plot(t, w_X2[,causalComp], type="h", col="darkblue", ann = FALSE, cex.axis=2)
title(xlab="SNP", cex.lab = 2)
title(ylab="SNP weight", cex.lab = 2)
lines(t[highsnp], w_X2[highsnp,causalComp], type="h", col="red", lwd=1.5)
points(t, w_X2[,causalComp], pch=16, col="darkblue")
points(t[highsnp], w_X2[highsnp,causalComp], pch=16, col="red")
legend("bottomright", X2names, col = c("red"), cex = 1, lty = 1, lwd=1.5)
title(main="SNP weight profile of PLSC", cex.main = 2)

dev.off()



# ----------------------------------------------
# Comparison of PLSC weights and PLSC-RP weights
# ----------------------------------------------

# load PLSC weights
load(paste(filenameRes, "w_X1_PLSC.RData", sep=""))
load(paste(filenameRes, "w_X2_PLSC.RData", sep=""))
w_X1_PLSC <- w_X1
w_X2_PLSC <- w_X2
w_X1_PLSC <- w_X1_PLSC[,causalComp]
w_X2_PLSC <- w_X2_PLSC[,causalComp]

# load PLSC-RP weights
load(paste(filenameRes, "w_X1_PLSC_RP.RData", sep=""))
load(paste(filenameRes, "w_X2_PLSC_RP.RData", sep=""))
w_X1_PLSC_RP <- w_X1
w_X2_PLSC_RP <- w_X2
w_X1_PLSC_RP <- w_X1_PLSC_RP[,causalComp]
w_X2_PLSC_RP <- w_X2_PLSC_RP[,causalComp]

# Similarity for voxel weights
SimX1 <- similarity(w_X1_PLSC, w_X1_PLSC_RP)

print("Similarity for X1")
print("=================")
print(paste("Pearson Correlation:", SimX1$sim_corr))
print(paste("Cosine Measure:", SimX1$sim_cos))
print(paste("extended Jaccard Similarity:", SimX1$sim_eJac))

# Similarity for SNP weights
SimX2 <- similarity(w_X2_PLSC, w_X2_PLSC_RP)

print("Similarity for X2")
print("=================")
print(paste("Pearson Correlation:", SimX2$sim_corr))
print(paste("Cosine Measure:", SimX2$sim_cos))
print(paste("extended Jaccard Similarity:", SimX2$sim_eJac))
