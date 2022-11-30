###############################
##
## Project: TreeKernel
##
## Purpose: Main TreeKernel function for data analysts to use
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2022-11-30
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

source("~/Documents/Research/Current/MetaboGuru/Carpenter/RCode/TreeKernels/gitFiles/ScoreSimFunctions.R")

## Function to perform TreeKernel test

## Y == The outcome of interest
## clinicalDat == The relevant clinical variables that will be used for clustering
## omicsData == The omics data used in the kernel association test
## graph == The graph used rescale omicsData if applicable. If left NULL, no rescaling is done
## nb.clust == The number of clusters expected from HCPC. See HCPC in the FactoMineR package for details
## classifier == The method used for embedding. If you are clustering on only continuous variables,
##                use 'continuous' principle components. Otherwise use 'discrete' for factor analysis
## ncp == The number of components taken from embedding and used for clustering
## outcomeType == "C" for continuous and "D" for dichotomous
## ... == Extra arguments passed to the HCPC function in FactoMineR

TreeKernel <- function(Y, clinicalData,  omicsData, graph=NULL,
                       nb.clust=0, classifier='continuous', 
                       ncp=5, outcomeType=c("C", "D"),
                       ...){
  
  covs <- clinicalData; Z <- omicsData
  dd <- data.frame(Y, covs)
  
  if(!is.null(graph)){
    I <- diag(length(V(graph)))
    NL <- as.matrix( graph.laplacian(graph, normalized = T) )
    RNL <- solve(I + tau*NL)
  } else{
    RNL <- diag(ncol(Z))
  }
  
  
  if(classifier=='continuous'){
    pca <- PCA(covs, graph=F, ncp=ncp)
    hcpc <- HCPC(pca, nb.clust=nb.clust, graph=F, ...)
  }
  if(classifier=='discrete'){
    famd <- FAMD(covs, graph=F, ncp=ncp)
    hcpc <- HCPC(famd, nb.clust=nb.clust, graph=F, ...)
  }
  
  dd <- data.frame(dd, hcpc$data.clust)
  
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
    
    if(outcomeType=="C"){
      pval[i] <- SKAT.c(Y~1, data=ds, K=Ks)$Qq
    }
    
    if(outcomeType=="D"){
      pval[i] <- SKAT.b(Y~1, data=ds, K=Ks)$p.value
    }
    
  }
  
  
  list(pvalHC = pval, tr = hcpc$call$t$tree)
}









