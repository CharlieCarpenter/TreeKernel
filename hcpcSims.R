###############################
##
## Project: MetaboGuru
##
## Purpose: Messing with hierarchical Simulations
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2022-05-23
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

## Helpful functions
source('~/Documents/Research/Current/MetaboGuru/Carpenter/RCode/ScoreSimFunctions.R')

# HCTree ------------------------------------------------------------------

## Function to simulate higher order interactions
## between covariate and kernel space
## From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6415770/
## graph for Z, Covariates, depth of kdTree, a numerical vector (length=2),
## sample size, size of graph, covariate beta vector,
## standard deviation of y
hcTreeIntDat <- function(graph, covs, zz, groups=2,
                         n, p, b0, sd.y, classifier='continuous'){
  
  stopifnot("All length(zz) must equal p"= all(sapply(zz,length)==p))
  stopifnot("length(zz) must equal groups"= length(zz)==groups)
  
  if(classifier=='continuous'){
    pca <- PCA(covs, graph=F)
    hcpcFull <- HCPC(pca, nb.clust=groups, min=groups, graph=F)
  }
  if(classifier=='discrete'){
    famd <- FAMD(covs, graph=F)
    hcpcFull <- HCPC(famd, nb.clust=4, min=4, graph=F)
  }
  
  tr <- hcpcFull$call$t$tree
  hcpc <- hcpcFull$data.clust
  hcpc$clust <- as.numeric(hcpc$clust)
  
  A <- as.matrix(get.adjacency(graph))
  test <- A + diag(p)
  Omega1 <- Danaher_pos_def(test)
  Z <- mvrnorm(n=nrow(covs), rep(0,p), solve(Omega1))
  
  ## Partitions
  if(any(hcpc$clust>groups)){
    hcpc$clust <- ifelse(hcpc$clust>groups,groups,hcpc$clust)
  }

  mX <- model.matrix(~X1+X2+X3+X4, data=covs)
  null <- mX%*%rep(b0, ncol(mX))
  Y <- VY <- ZK <- numeric(n)
  
  for(i in 1:length(unique(hcpc$clust)) ){
    kk <- hcpc$clust==i 
    Zk <- Z[kk,] %*% zz[[i]]
    ZK[kk] <- Zk
    Y[kk] <- null[kk] + Zk
  }
  
  Ve <- rnorm(n=n, sd=sd.y)
  Y <- Y + Ve
  
  list(Y=Y, graph=graph, covs=covs, Z=Z, tr=tr, Ve=var(Ve), VY=var(ZK))
}


## Function to perform TreeKernel on output from hcTreeIntDat
## This is the function used for simulations
hcTreeTest <- function(hcDat, formula.H0, groups=2, classifier='continuous', ncp=5){
  
  Y <- hcDat$Y; covs <- hcDat$covs
  graph <- hcDat$graph; Z <- hcDat$Z
  dd <- data.frame(Y, covs)
  mod <- lm(formula.H0, data=dd)
  R <- mod$residuals
  dd <- data.frame(R,dd)
  
  I <- diag(length(V(graph)))
  NL <- as.matrix( graph.laplacian(graph, normalized = T) )
  RNL <- solve(I + tau*NL)
  
  if(classifier=='continuous'){
    pca <- PCA(covs, graph=F, ncp=ncp)
    hcpc <- HCPC(pca, nb.clust = groups, min=groups, graph=F)$data.clust
  }
  if(classifier=='discrete'){
    famd <- FAMD(covs, graph=F, ncp=ncp)
    hcpc <- HCPC(famd, nb.clust = groups, min=groups, graph=F)$data.clust
  }
  
  dd <- data.frame(dd, hcpc)
  
  ## testing within each partition
  pval <- numeric()
  for(i in 1:length(unique(hcpc$clust)) ){
    kk <- dd$clust==i
    if(sum(kk)==1) next
    ds <- dd[kk,]
    Zs <- scale(Z[kk, ,drop=F])
    rho <- median(dist(Zs))
    Zl <- try(Zs %*% RNL, silent = T)
    
    Ks <- Gaussian_kernel(rho, Zl)
    
    pval[i] <- SKAT.c(R~1, data=ds, K=Ks)$Qq
  }
  
  Zs <- scale(Z)
  rho <- median(dist(Zs))
  Zl <- Zs %*% RNL
  K <- Gaussian_kernel(rho, Zl)
  
  fullP <- SKAT.c(formula.H0, data=dd, K=K)$Qq
  
  list(pvalHC = pval, fullP=fullP, tr = hcDat$tr)
}




hcTreeIntDatOneTree <- function(graph, covs, hcpcFull, zz, groups=2,
                                n, p, b0, sd.y){
  
  stopifnot("All length(zz) must equal p"= all(sapply(zz,length)==p))
  stopifnot("length(zz) must equal groups"= length(zz)==groups)
  
  tr <- hcpcFull$call$t$tree
  hcpc <- hcpcFull$data.clust
  hcpc$clust <- as.numeric(hcpc$clust)
  
  A <- as.matrix(get.adjacency(graph))
  test <- A + diag(p)
  Omega1 <- Danaher_pos_def(test)
  Z <- mvrnorm(n=nrow(covs), rep(0,p), solve(Omega1))
  
  ## Partitions
  if(any(hcpc$clust>groups)){
    hcpc$clust <- ifelse(hcpc$clust>groups,groups,hcpc$clust)
  }
  
  mX <- as.matrix(cbind(1, covs))
  null <- mX%*%b0 
  Y <- VY <- ZK <- numeric(n)
  for(i in 1:groups){
    kk <- hcpc$clust==i 
    Zk <- Z[kk,] %*% zz[[i]]
    
    ZK[kk] <- Zk
    Y[kk] <- null[kk] + Zk
  }
  
  Ve <- rnorm(n=n, sd=sd.y)
  Y <- Y + Ve
  
  list(Y=Y, graph=graph, covs=covs, Z=Z, tr=tr, Ve=var(Ve), VY=var(ZK))
}




