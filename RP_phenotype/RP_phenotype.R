# _______________________________________________________________________
# Using PLSC-RP for optimized analysis of high-dimensional phenotype data
# _______________________________________________________________________


# Clear the environment
rm(list=ls())
library(far)
library(proxy)


###################################################################################################################################
# insert filename of example data here


# Example Data with 1000 voxels

# fMRI data (phenotype data)
filenameX1 <- "C:/Users/cgrellmann/Documents/Privat/Promotion/Paper/Frontiers/PlscRP/RP_phenotype/ExampleData_1000voxel/X1.RData"
# SNP data (genotype data)
filenameX2 <- "C:/Users/cgrellmann/Documents/Privat/Promotion/Paper/Frontiers/PlscRP/RP_phenotype/ExampleData_1000voxel/X2.RData"

# ID of causal voxel
filenameID_causalX1 <- "C:/Users/cgrellmann/Documents/Privat/Promotion/Paper/Frontiers/PlscRP/RP_phenotype/ExampleData_1000voxel/IDcausalX1.RData"
# ID of voxels collinear to the causal voxel (r>=0.8)
filenameID_collinearX1 <- "C:/Users/cgrellmann/Documents/Privat/Promotion/Paper/Frontiers/PlscRP/RP_phenotype/ExampleData_1000voxel/IDcollinearX1.RData"
# ID of causal SNPs
filenameID_causalX2 <- "C:/Users/cgrellmann/Documents/Privat/Promotion/Paper/Frontiers/PlscRP/RP_phenotype/ExampleData_1000voxel/IDcausalX2.RData"


# insert filename of results here
filenameRes <- "C:/Users/cgrellmann/Documents/Privat/Promotion/Paper/Frontiers/PlscRP/RP_phenotype/ExampleData_1000voxel/"


# orthogonalization of random matrix using Gram-Schmidt algorithm
orthonorm <- TRUE


# number of permutations
nperms <- 5000
# more accurate results for high number of permutations
# faster results for small number of permutations


# set working directory
setwd("C:/Users/cgrellmann/Documents/Privat/Promotion/Paper/Frontiers/PlscRP/functions")

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

load(filenameX1)               # fMRI data (phenotype data)
load(filenameX2)               # SNP data (genotype data)
load(filenameID_causalX1)      # ID of causal voxel
load(filenameID_collinearX1)   # ID of voxels collinear to the causal voxel (r>=0.8)
load(filenameID_causalX2)      # ID of causal SNPs



#-------------------------------------
# Check dimensionality of example data
#-------------------------------------

N <- nrow(X1)
D1 <- ncol(X1)
D2 <- ncol(X2)

print("Dimensionality of example data")
print("==============================")
print(paste("Sample size:", N))
print(paste("Number of voxels:", D1))
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
# RP for dimensionality reduction in X1
# -------------------------------------

# define dimensionality of low-dimensional subspace using Menon theorem
epsi <- 1.0    # distance preservation factor
prob <- 0.05   # 1 - probability of distance preservation

kMen <- lbMenon(N, epsi, prob)
kMen <- kMen$k

if(ceiling(kMen) < N){
  K <- N
}else{
  K <- ceiling(kMen)
}

print("Menon theorem")
print("=============")
print(paste("Menon lower bound:", kMen))
print(paste("Specified dimensionality of low-dimensional subspace:", K))


# create i.i.d. N(0,1) matrix
randmat <- matrix(rnorm(D1*K, mean=0, sd=1), D1, K)   # rows and columns are N(0,1)-distributed

# orthogonalization using Gram-Schmidt algorith
if(orthonorm){
  randmat <- orthonormalization(randmat, basis=FALSE, norm=TRUE)
}

# dimension reduction using RP
X1.RP <- X1.zscore %*% randmat



#------------
# Run PLSC-RP
#------------

# Permutation test
comp <- 1        # component for which permutation test should be performed

perm.PLSC_RP <- PLSC.permute(X1.RP, X2.zscore, nperms, comp)

print("Permutation test for PLSC-RP")
print("============================")
print(paste("Covariance of latent variables:", perm.PLSC_RP$cov))
print(paste("p-value of covariance:", perm.PLSC_RP$pval_cov))


# PLSC-RP analysis
out.PLSC_RP <- PLSC(X1.RP, X2.zscore)

w_X1.RP <- out.PLSC_RP$A
w_X2 <- out.PLSC_RP$B
LVX1.RP <- out.PLSC_RP$LV_X
LVX2 <- out.PLSC_RP$LV_Y
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

w_X1.RP <- w_X1.RP[,1:sel]
w_X2 <- w_X2[,1:sel]
LVX1.RP <- LVX1.RP[,1:sel]
LVX2 <- LVX2[,1:sel]
covval <- covval[1:sel]


print("Results of PLSC-RP")
print("==================")
print(paste("Number of components:", sel))
print(paste("Variance explained:", sum(totalVar[1:sel])))
print(paste("Covariance of latent variables:", covval))



#------------------------------------------
# Out-of-sample covariance using 10-fold CV
#------------------------------------------

cov_oos_PLSC_RP <- Cov_oos_CV(X1.RP, X2.zscore, sel)

cov_oos_mean <- cov_oos_PLSC_RP$cov_oos_mean
causalComp <- which.max(diff(cov_oos_mean))
raw_cov_oos <- diff(cov_oos_mean)

print("Out-of-sample covariance using 10-fold CV")
print("=========================================")
print(paste("Causal Component:", causalComp))
print(paste("Out-of-sample covariance of causal component:", raw_cov_oos[causalComp]))



#-------------------------------------------------
# Backtransformation to get original voxel weights
#-------------------------------------------------

covval <- covval[1]
w_X2 <- w_X2[,1]
w_X2 <- matrix(w_X2)
LVX2 <- LVX2[,1]
LVX2 <- matrix(LVX2)

w_X1 <- (1/(covval * sum(w_X2^2))) * rX1X2 %*% w_X2
colnorm <- sqrt(sum(w_X1^2))
w_X1 <- (1/colnorm) * w_X1
LVX1 <- LVX1.RP[,1]
LVX1 <- matrix(LVX1)



#--------------------------
# Saving results of PLSC-RP
#--------------------------
numVox <- nrow(w_X1)

save(w_X1, file = paste(filenameRes, "w_X1_PLSC_RP_", numVox, "voxel.RData", sep=""))
save(w_X2, file = paste(filenameRes, "w_X2_PLSC_RP_", numVox, "voxel.RData", sep=""))
save(LVX1, file = paste(filenameRes, "LVX1_PLSC_RP_", numVox, "voxel.RData", sep=""))
save(LVX2, file = paste(filenameRes, "LVX2_PLSC_RP_", numVox, "voxel.RData", sep=""))



#--------------------------
# Presenting results for X1
#--------------------------

t <- 1:nrow(w_X1)

# entire voxel weight profile
graphics.off()
X11()

dev.copy(pdf, file = paste(filenameRes, "w_X1_PLSC_RP_", numVox, "voxel_total.pdf", sep=""))

par(mar=c(5.1,4.5,4.1,2.1))
plot(t, w_X1, type="h", col="darkblue", ann = FALSE, cex.axis=2)
title(xlab="voxel", cex.lab = 2)
title(ylab="voxel weight", cex.lab = 2)
lines(t[IDcollinearX1], w_X1[IDcollinearX1], type="h", col="gold", lwd=1.5)
lines(t[IDcausalX1], w_X1[IDcausalX1], type="h", col="red", lwd=1.5)
legend("bottomright", c("causal voxel","causal voxels r=0.8"), col = c("red","gold"), cex = 1, lty = 1, lwd=c(1.5))
title(main="Voxel weight profile of PLSC-RP", cex.main = 2)

dev.off()

# zoomed in profile presented in the manuscript
graphics.off()
X11()

dev.copy(pdf, file = paste(filenameRes, "w_X1_PLSC_RP_", numVox, "voxel.pdf", sep=""))

par(mar=c(5.1,4.5,4.1,2.1))
plot(t[(IDcausalX1-30):(IDcausalX1+70)], w_X1[(IDcausalX1-30):(IDcausalX1+70)], type="h", col="darkblue", ann = FALSE, cex.axis=2)
title(xlab="voxel", cex.lab = 2)
title(ylab="voxel weight", cex.lab = 2)
lines(t[IDcollinearX1], w_X1[IDcollinearX1], type="h", col="gold", lwd=1.5)
lines(t[IDcausalX1], w_X1[IDcausalX1], type="h", col="red", lwd=1.5)
legend("bottomright", c("causal voxel","causal voxels r=0.8"), col = c("red","gold"), cex = 1, lty = 1, lwd=c(1.5))
title(main="Voxel weight profile of PLSC-RP", cex.main = 2)

dev.off()



#--------------------------
# Presenting results for X2
#--------------------------

t <- 1:nrow(w_X2)

graphics.off()
X11()

dev.copy(pdf, file = paste(filenameRes, "w_X2_PLSC_RP_", numVox, "voxel.pdf", sep=""))

par(mar=c(5.1,4.5,4.1,2.1))
plot(t, w_X2, type="h", col="darkblue", ann = FALSE, cex.axis=2)
title(xlab="SNP", cex.lab = 2)
title(ylab="SNP weight", cex.lab = 2)
lines(t[IDcausalX2], w_X2[IDcausalX2], type="h", col="red", lwd=1.5)
points(t, w_X2, pch=16, col="darkblue")
points(t[IDcausalX2], w_X2[IDcausalX2], pch=16, col="red")
legend("bottomleft", c("causal SNPs"), col = c("red"), cex = 1, lty = 1, lwd=1.5)
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

w_X1 <- w_X1[,1]
w_X2 <- w_X2[,1]
LVX1 <- LVX1[,1]
LVX2 <- LVX2[,1]

w_X1 <- matrix(w_X1)
w_X2 <- matrix(w_X2)
LVX1 <- matrix(LVX1)
LVX2 <- matrix(LVX2)

covval <- covval[1]



#--------------------------
# Saving results of PLSC-RP
#--------------------------

save(w_X1, file = paste(filenameRes, "w_X1_PLSC_", numVox, "voxel.RData", sep=""))
save(w_X2, file = paste(filenameRes, "w_X2_PLSC_", numVox, "voxel.RData", sep=""))
save(LVX1, file = paste(filenameRes, "LVX1_PLSC_", numVox, "voxel.RData", sep=""))
save(LVX2, file = paste(filenameRes, "LVX2_PLSC_", numVox, "voxel.RData", sep=""))



#--------------------------
# Presenting results for X1
#--------------------------

t <- 1:nrow(w_X1)

# entire voxel weight profile
graphics.off()
X11()

dev.copy(pdf, file = paste(filenameRes, "w_X1_PLSC_", numVox, "voxel_total.pdf", sep=""))

par(mar=c(5.1,4.5,4.1,2.1))
plot(t, w_X1, type="h", col="darkblue", ann = FALSE, cex.axis=2)
title(xlab="voxel", cex.lab = 2)
title(ylab="voxel weight", cex.lab = 2)
lines(t[IDcollinearX1], w_X1[IDcollinearX1], type="h", col="gold", lwd=1.5)
lines(t[IDcausalX1], w_X1[IDcausalX1], type="h", col="red", lwd=1.5)
legend("bottomright", c("causal voxel","causal voxels r=0.8"), col = c("red","gold"), cex = 1, lty = 1, lwd=c(1.5))
title(main="Voxel weight profile of PLSC", cex.main = 2)

dev.off()

# zoomed in profile presented in the manuscript
graphics.off()
X11()

dev.copy(pdf, file = paste(filenameRes, "w_X1_PLSC_", numVox, "voxel.pdf", sep=""))

par(mar=c(5.1,4.5,4.1,2.1))
plot(t[(IDcausalX1-30):(IDcausalX1+70)], w_X1[(IDcausalX1-30):(IDcausalX1+70)], type="h", col="darkblue", ann = FALSE, cex.axis=2)
title(xlab="voxel", cex.lab = 2)
title(ylab="voxel weight", cex.lab = 2)
lines(t[IDcollinearX1], w_X1[IDcollinearX1], type="h", col="gold", lwd=1.5)
lines(t[IDcausalX1], w_X1[IDcausalX1], type="h", col="red", lwd=1.5)
legend("bottomright", c("causal voxel","causal voxels r=0.8"), col = c("red","gold"), cex = 1, lty = 1, lwd=c(1.5))
title(main="Voxel weight profile of PLSC", cex.main = 2)

dev.off()



#--------------------------
# Presenting results for X2
#--------------------------

t <- 1:nrow(w_X2)

graphics.off()
X11()

dev.copy(pdf, file = paste(filenameRes, "w_X2_PLSC_", numVox, "voxel.pdf", sep=""))

par(mar=c(5.1,4.5,4.1,2.1))
plot(t, w_X2, type="h", col="darkblue", ann = FALSE, cex.axis=2)
title(xlab="SNP", cex.lab = 2)
title(ylab="SNP weight", cex.lab = 2)
lines(t[IDcausalX2], w_X2[IDcausalX2], type="h", col="red", lwd=1.5)
points(t, w_X2, pch=16, col="darkblue")
points(t[IDcausalX2], w_X2[IDcausalX2], pch=16, col="red")
legend("bottomleft", c("causal SNPs"), col = c("red"), cex = 1, lty = 1, lwd=1.5)
title(main="SNP weight profile of PLSC", cex.main = 2)

dev.off()



# ----------------------------------------------
# Comparison of PLSC weights and PLSC-RP weights
# ----------------------------------------------

# load PLSC weights
load(paste(filenameRes, "w_X1_PLSC_", numVox, "voxel.RData", sep=""))
load(paste(filenameRes, "w_X2_PLSC_", numVox, "voxel.RData", sep=""))
w_X1_PLSC <- w_X1
w_X2_PLSC <- w_X2

# load PLSC-RP weights
load(paste(filenameRes, "w_X1_PLSC_RP_", numVox, "voxel.RData", sep=""))
load(paste(filenameRes, "w_X2_PLSC_RP_", numVox, "voxel.RData", sep=""))
w_X1_PLSC_RP <- w_X1
w_X2_PLSC_RP <- w_X2

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
