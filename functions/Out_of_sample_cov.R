# _________________________________________
# Out-of-sample covariance using 10-fold CV
# _________________________________________

# input:
# ------
# X1 = N x D1 phenotype matrix
# X2 = N x D2 genotype matrix
# sel = number of components that explain 80% of variance

# output:
# -------
# cov_oos = cummulative out-of-sample covariance for the 10 folds
# cov_oos_mean = mean out-of-sample covariance



Cov_oos_CV <- function(X1, X2, sel){
  
  N <- nrow(X1)
  shuffle <- sample(N)
  
  X1_s <- X1[shuffle,]
  X2_s <- X2[shuffle,]
  
  fold <- floor(N/10)
  
  f1 <- 1:fold
  f2 <- (fold+1):(2*fold)
  f3 <- (2*fold+1):(3*fold)
  f4 <- (3*fold+1):(4*fold)
  f5 <- (4*fold+1):(5*fold)
  f6 <- (5*fold+1):(6*fold)
  f7 <- (6*fold+1):(7*fold)
  f8 <- (7*fold+1):(8*fold)
  f9 <- (8*fold+1):(9*fold)
  f10 <- (9*fold+1):N
  
  X1train <- vector("list", length=10)
  X1test <- vector("list", length=10)
  X2train <- vector("list", length=10)
  X2test <- vector("list", length=10)
  cov_oos_CV <- vector("list", length=10)
  
  X1train[[1]] <- X1_s[c(f2,f3,f4,f5,f6,f7,f8,f9,f10),]
  X1test[[1]] <- X1_s[f1,]
  X2train[[1]] <- X2_s[c(f2,f3,f4,f5,f6,f7,f8,f9,f10),]
  X2test[[1]] <- X2_s[f1,]
  
  X1train[[2]] <- X1_s[c(f1,f3,f4,f5,f6,f7,f8,f9,f10),]
  X1test[[2]] <- X1_s[f2,]
  X2train[[2]] <- X2_s[c(f1,f3,f4,f5,f6,f7,f8,f9,f10),]
  X2test[[2]] <- X2_s[f2,]
  
  X1train[[3]] <- X1_s[c(f1,f2,f4,f5,f6,f7,f8,f9,f10),]
  X1test[[3]] <- X1_s[f3,]
  X2train[[3]] <- X2_s[c(f1,f2,f4,f5,f6,f7,f8,f9,f10),]
  X2test[[3]] <- X2_s[f3,]
  
  X1train[[4]] <- X1_s[c(f1,f2,f3,f5,f6,f7,f8,f9,f10),]
  X1test[[4]] <- X1_s[f4,]
  X2train[[4]] <- X2_s[c(f1,f2,f3,f5,f6,f7,f8,f9,f10),]
  X2test[[4]] <- X2_s[f4,]
  
  X1train[[5]] <- X1_s[c(f1,f2,f3,f4,f6,f7,f8,f9,f10),]
  X1test[[5]] <- X1_s[f5,]
  X2train[[5]] <- X2_s[c(f1,f2,f3,f4,f6,f7,f8,f9,f10),]
  X2test[[5]] <- X2_s[f5,]
  
  X1train[[6]] <- X1_s[c(f1,f2,f3,f4,f5,f7,f8,f9,f10),]
  X1test[[6]] <- X1_s[f6,]
  X2train[[6]] <- X2_s[c(f1,f2,f3,f4,f5,f7,f8,f9,f10),]
  X2test[[6]] <- X2_s[f6,]
  
  X1train[[7]] <- X1_s[c(f1,f2,f3,f4,f5,f6,f8,f9,f10),]
  X1test[[7]] <- X1_s[f7,]
  X2train[[7]] <- X2_s[c(f1,f2,f3,f4,f5,f6,f8,f9,f10),]
  X2test[[7]] <- X2_s[f7,]
  
  X1train[[8]] <- X1_s[c(f1,f2,f3,f4,f5,f6,f7,f9,f10),]
  X1test[[8]] <- X1_s[f8,]
  X2train[[8]] <- X2_s[c(f1,f2,f3,f4,f5,f6,f7,f9,f10),]
  X2test[[8]] <- X2_s[f8,]
  
  X1train[[9]] <- X1_s[c(f1,f2,f3,f4,f5,f6,f7,f8,f10),]
  X1test[[9]] <- X1_s[f9,]
  X2train[[9]] <- X2_s[c(f1,f2,f3,f4,f5,f6,f7,f8,f10),]
  X2test[[9]] <- X2_s[f9,]
  
  X1train[[10]] <- X1_s[c(f1,f2,f3,f4,f5,f6,f7,f8,f9),]
  X1test[[10]] <- X1_s[f10,]
  X2train[[10]] <- X2_s[c(f1,f2,f3,f4,f5,f6,f7,f8,f9),]
  X2test[[10]] <- X2_s[f10,]
  
  for(i in 1:10){
    X1tr <- X1train[[i]]
    X2tr <- X2train[[i]]
    X1te <- X1test[[i]]
    X2te <- X2test[[i]]
    
    cat("PLSC for CV fold", i, "\n")
    
    PLSC_CV <- PLSC(X1tr, X2tr)
    
    w_X1 <- PLSC_CV$A
    w_X2 <- PLSC_CV$B
    
    Z_X1 <- X1te %*% w_X1
    Z_X2 <- X2te %*% w_X2
    
    # out of sample covariance
    cov_oos <- matrix(0,1,ncol(w_X2))
    
    cov <- cov(Z_X1,Z_X2)
    cov <- diag(cov)
    
    for(j in 1:ncol(w_X2)){
      cov_oos[j] <- sum(cov[1:j])
    }
    
    cov_oos_CV[[i]] <- cov_oos
  }
  
  # mean of cov_oos of all CV
  cov_oos <- rbind(cov_oos_CV[[1]],cov_oos_CV[[2]],cov_oos_CV[[3]],cov_oos_CV[[4]],cov_oos_CV[[5]],cov_oos_CV[[6]],cov_oos_CV[[7]],cov_oos_CV[[8]],cov_oos_CV[[9]],cov_oos_CV[[10]])
  cov_oos_mean <- colMeans(cov_oos)
  
  # selected columns explaining 80% of variance
  cov_oos <- cov_oos[,1:sel]
  cov_oos_mean <- cov_oos_mean[1:sel]
  
  # include zero
  cov_oos <- cbind(matrix(c(0,0,0,0,0,0,0,0,0,0), ncol=1), cov_oos)
  cov_oos_mean <- c(0, cov_oos_mean)
  
  return(list(cov_oos=cov_oos, cov_oos_mean=cov_oos_mean))
}
