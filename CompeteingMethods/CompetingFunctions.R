###############################
##
## Project: MetaboGuru
##
## Purpose: TreeKernel Competitor simulations
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2022-08-03
##
## ---------------------------
## Notes:
##
## ---------------------------

# Helpful functions
`%nin%` <- Negate(`%in%`)
library(tidyverse); library(magrittr)

## Data Read In ----

## Test for relationships within the covariate partitions
## kdDat is the output from "kdTreeIntDat" function
pcaTreeTest <- function(hcDat, formula.H0, groups=2, classifier='continuous', ncp=5){
  
  Y <- hcDat$Y; covs <- hcDat$covs
  graph <- hcDat$graph; Z <- hcDat$Z
  dd <- data.frame(Y, covs)
  mod <- lm(formula.H0, data=dd)
  
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
  dd <- data.frame(Y=Y, hcpc)
  
  ## testing within each partition
  pval <- numeric()
  for(i in 1:length(unique(hcpc$clust)) ){
    kk <- dd$clust==i
    if(sum(kk)==1) next
    ds <- dd[kk,]
    Zs <- scale(Z[kk, ,drop=F])
    
    pval[i] <- CCpcaTest(formula.H0, ds, Zs, m=5)
  }
  
  Zs <- scale(Z)
  rho <- median(dist(Zs))
  Zl <- Zs %*% RNL
  K <- Gaussian_kernel(rho, Zl)
  
  fullP <- SKAT.c(formula.H0, data=dd, K=K)$Qq
  
  list(pvalHC = pval, fullP=fullP, tr = hcDat$tr)
}


## Test for relationships within the covariate partitions
## kdDat is the output from "kdTreeIntDat" function
simesTreeTest <- function(hcDat, formula.H0, groups=2, classifier='continuous', ncp=5){
  
  Y <- hcDat$Y; covs <- hcDat$covs
  graph <- hcDat$graph; Z <- hcDat$Z
  dd <- data.frame(Y, covs)
  mod <- lm(formula.H0, data=dd)
  
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
  dd <- data.frame(Y=Y, hcpc)
  
  ## testing within each partition
  pval <- numeric()
  for(i in 1:length(unique(hcpc$clust)) ){
    kk <- dd$clust==i
    if(sum(kk)==1) next
    ds <- dd[kk,]
    Zs <- scale(Z[kk, ,drop=F]) 
    Zl <- try(Zs %*% RNL, silent = T)
    
    pval[i] <- Simes(formula.H0, ds, .Z=Zl)
  }
  
  Zs <- scale(Z)
  Zl <- Zs %*% RNL
  
  fullP <- Simes(formula.H0, dd, .Z=Zl)
  
  list(pvalHC = pval, fullP=fullP, tr = hcDat$tr)
}


# Core Functions ----------------------------------------------------------

CCpcaTest <- function(formula.H0, ds, Zs, m){
  SVD <- svd(Zs)
  scores <- SVD$u %*% diag(SVD$d)
  scrs <- scores[,1:m]
  
  ll1 <- lm(formula.H0, data=ds)
  dss <- cbind(ds, scrs)
  
  fz <- update(formula.H0, ~ . + `1`+`2`+`3`+`4`+`5`)
  
  ll2 <- lm(fz, data=dss)
  
  anova(ll1, ll2)$`Pr(>F)`[2]
}

Simes <- function(formula.H0, dd, .Z){
  m <- ncol(.Z)
  ps <- numeric(m)
  
  formula.Z <- update(formula.H0,  ~ . + Z)
  
  ## Univariate tests
  for(i in 1:m){
    dz <- data.frame(dd, Z = .Z[,i])
    ll <- lm(formula.Z, data = dz)
    ps[i] <- coef(summary(ll))[-1, 4]
  }
  
  min(sort(ps)*(m)/(1:m)) # simes p-value adjustment
}
