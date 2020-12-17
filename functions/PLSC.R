# 
# Plsc-RP: Partial Least Squares Correlation - Random Projection: a method to efficiently solve high-dimensional multimodal problems 
# Copyright (C) 2016. Claudia Grellmann & The O'BRAIN Lab
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# The code was written by Claudia Grellmann. If you have questions on the method and/ or data, please get in touch with The O'BRAIN Lab via its website: https://obrainlab.com/
#
#
#_________________________________
# Partial Least Squares Correlation
# _________________________________

# ----------------
# Permutation test
# ----------------

# input:
# ------
# X = N x m phenotype matrix
# Y = N x n genotype matrix
# nperms = number of permutations
# comp = component for which permutation test should be performed

# output:
# -------
# cov = covariance of original scores
# covperms = covariance of permuted scores
# pval_cov = p value of covariance
# nperms = number of permutations

PLSC.permute <- function(X, Y, nperms, comp){
    if(ncol(X) < 2) stop("Need at least 2 features in data set X.")
    if(ncol(Y) < 2) stop("Need at least 2 features in data set Y.")
    if(sum(is.na(X)) + sum(is.na(Y)) > 0) stop("Cannot have NAs in X or Y")
    if(nrow(X) != nrow(Y)) stop("X and Y must have same number of rows")
    
    out <- PLSC(X, Y)
    d <- out$d[comp]   # covariance of scores
    
    dperms <- apply(matrix(1:nperms), 1, function(n){
      return(PLSC(X, Y[sample(1:nrow(Y)),])$d[comp])
    })
    
    # p value of covariance of scores
    d <- as.numeric(d)
    dperms <- as.numeric(dperms)
    
    pval_cov <- mean(dperms >= d)
    
    results <- list(cov=d, covperms=dperms, pval_cov=pval_cov, nperms=nperms)
    return(results)
    
}



# -------------
# PLSC analysis
# -------------

# input:
# ------
# X = N x m phenotype matrix
# Y = N x n genotype matrix

# output:
# -------
# A = left singular vector
# B = right singular vector
# d = covariance between latent variables
# LV_X = latent variables of X
# LV_Y = latent variables of Y
# p - min(m,n)
# totalVar - proportion of the sum of squared cross-block correlations explained by A and B

PLSC <- function(X, Y){
    if(ncol(X) < 2) stop("Need at least 2 features in data set X.")
    if(ncol(Y) < 2) stop("Need at least 2 features in data set Y.")
    if(sum(is.na(X)) + sum(is.na(Y)) > 0) stop("Cannot have NAs in X or Y")
    if(nrow(X) != nrow(Y)) stop("X and Y must have same number of rows")
  
    R_XY <- cor(X, Y)   # cross correlation matix m x n
    
    s <- svd(R_XY)      # singular value decomposition of cross correlation matix R_XY

    d <- s$d            # covariance between latent variables (length min(m,n))
    A <- s$u            # left singular vector (dimension m x min(m,n))
    B <- s$v            # right singular vector (dimension n x min(m,n))

    p <- length(d)
    
    # computation of latent variables
    LV_X <- X %*% A     # dimension N x min(m,n)
    LV_Y <- Y %*% B     # dimension N x min(m,n)

    # proportion of the sum of squared cross-block correlations explained by A and B
    totalVar <- numeric(p)
    for(i in 1:p){
      totalVar[i] <- d[i]^2 / sum(d^2)
    }
    totalVar <- totalVar * 100

    return(list(A=A, B=B, d=d, LV_X=LV_X, LV_Y=LV_Y, p=p, totalVar=totalVar))
}