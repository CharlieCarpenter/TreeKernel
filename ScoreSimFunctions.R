###############################
##
## Project: MetaboGuru
##
## Purpose: Score Stat Simulation functions
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2020-11-10
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

# Helpful functions
# Helpful functions
`%nin%` <- Negate(`%in%`)

pkgs <- c("tidyverse", "magrittr", "igraph", "matrixcalc",
          "MASS", "diffusr", "CompQuadForm", "FactoMineR")
suppressMessages(lapply(pkgs, library, character.only = T))

# Functions ---------------------------------------------------------------

## Function to make matrix positive definite following approach of Danaher et al (2014)
Danaher_pos_def <- function(m, cc = 3.4, dd = 0.95){
  AA <- m*cc + t(m*cc)
  AA[AA>0] <- 1
  AA <- AA-diag(diag(AA))+diag(nrow(m))*dd
  AA <- AA/( as.matrix(rep(1,nrow(m))) ) %*% t(1.4*rowSums(abs(AA)))
  AA <- (AA+t(AA))/2
  AA <- AA-diag(diag(AA))+diag(nrow(m))*dd
  AA
}

## Trace of matrix
tr <- function(x){
  if(!is.matrix(x)) stop('x must be a matrix to find the trace')
  sum(diag(x))
}

Gaussian_kernel <- function(rho, Z){
  exp(-(1/rho)*as.matrix(dist(Z, method = "euclidean", upper = T)^2))
}

## Linear and poly Kern
plyKern <- function(Z, pow, rho=0){
  (Z%*%t(Z) + rho)^pow
}

## Data Generation
perfect_net_data <- function(graph, n, p, b0, sd.y, X, zz){
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A + diag(p)
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=n, rep(0,p), solve(Omega1))
    mX <- model.matrix(~X1+X2, data = X)
    
    Y <- rnorm(n, mX%*%b0+Z%*%zz, sd.y)
    # V.Y <- var(b0+Z%*%zz); V.e <- var(Y - b0+Z%*%zz)
    
    # total_edge_n <- sum(A==1)
  }
  list(Y=Y, Z=Z)
}

# Pieces for davies -------------------------------------------------

#Compute the tail probability of 1-DF chi-square mixtures
KAT.pval <- function(Q.all, lambda, acc=1e-9,lim=1e6){
  pval = rep(0, length(Q.all))
  i1 = which(is.finite(Q.all))
  for(i in i1){
    tmp <- davies(Q.all[i],lambda,acc=acc,lim=lim)
    pval[i] = tmp$Qq
    
    if(tmp$ifault>0) warning(paste("ifault =", tmp$ifault))
      # pval[i] = Sadd.pval(Q.all[i],lambda)
  }
  return(pval)
}

SKAT.c <- function(formula.H0, data = NULL, K, adjusted = T,
                   acc = 0.00001, lim = 10000, tol = 1e-10) {
  
  m0 <- lm(formula.H0, data)
  mX <- model.matrix(formula.H0, data)
  
  res <- resid(m0); df <- nrow(mX)-ncol(mX)
  s2 <- sum(res^2)/df
  
  P0  <- diag(nrow(mX)) - mX %*% (solve(t(mX) %*% mX) %*% t(mX))
  # print(dim(K));print(dim(P0))
  
  PKP <- P0 %*% K %*% P0
  
  if(adjusted){
    q <- as.numeric(res %*% K %*% res /(s2*df))
    ee <- eigen(PKP - q * P0, symmetric = T)
    q <- 0 ## Redefining for adjusted stat
  } else{
    q <- as.numeric(res %*% K %*% res / s2)
    ee <- eigen(PKP, symmetric = T) 
  }
  
  lambda <- ee$values[abs(ee$values) >= tol]
  dav <- davies(q, lambda = sort(lambda, decreasing=T),
                acc = acc, lim = lim)
  
  c(dav, Q.adj=q)
}

SKAT.b <- function(formula.H0, data = NULL, K, adjusted = F,
                   acc = 0.00001, lim = 10000, tol = 1e-10) {
  
  X1 <- model.matrix(formula.H0, data)
  lhs <- formula.H0[[2]]
  y <- eval(lhs, data)
  
  y <- factor(y)
  
  
  if (nlevels(y) != 2) {
    stop('The phenotype is not binary!\n')
  } else {
    y <- as.numeric(y) - 1
  }
  
  glmfit <- glm(y ~ X1 - 1, family = binomial)
  
  betas <- glmfit$coef
  mu  <- glmfit$fitted.values
  eta <- glmfit$linear.predictors
  res.wk <- glmfit$residuals
  res <- y - mu
  
  w   <- mu * (1-mu)
  sqrtw <- sqrt(w)
  
  adj <- sum((sqrtw * res.wk)^2) 
  
  DX12 <- sqrtw * X1
  
  qrX <- qr(DX12)
  Q <- qr.Q(qrX)
  Q <- Q[, 1:qrX$rank, drop=FALSE]
  
  P0 <- diag(nrow(X1)) - Q %*% t(Q)
  
  DKD <- tcrossprod(sqrtw) * K
  tQK <- t(Q) %*% DKD
  QtQK <- Q %*% tQK 
  PKP <- DKD - QtQK - t(QtQK) + Q %*% (tQK %*% Q) %*% t(Q)
  q <- as.numeric(res %*% K %*% res) / adj
  ee <- eigen(PKP - q * P0, symmetric = T, only.values=T)  		
  lambda <- ee$values[abs(ee$values) >= tol]
  
  p.value <- KAT.pval(0, lambda=sort(lambda, decreasing = T), acc = acc, lim = lim) 
  
  return(list(p.value=p.value, Q.adj = q))
}

# Useful functions --------------------------------------------------------

## Functions
count_shared_edges <- function(adjTrue1, adjTrue2, directed = T) {
  if(directed) se <- sum((adjTrue1 + adjTrue2)==2) / 2
  else se <- sum((adjTrue1 + adjTrue2)==2)
  se
}

## Function to randomly(ish) make a graph sparse
new.degs <- function(gg, p = 0.25){
  dd <- degree(gg)
  
  ## Replacing ~half the vertices with degree 1 with degree 0
  new.degs <- sapply(dd, function(deg, cut){
    if(deg <= cut) deg <- rbinom(1,1,0.5)
    deg
  }, cut = quantile(dd, probs = p)) %>%
    unlist
  
  new.degs
}

## Function to permute some percentage of edges
edge.perm <- function(adj, perc.perm){
  p <- nrow(adj) ## Number of vertices
  perm <- 1:p
  change <- sample(p, floor(p*perc.perm))
  
  change_idx <- change[order(change)]
  
  for (cidx in 1:length(change)){
    perm[change_idx[cidx]] <- change[cidx]
  }
  
  adj[perm, perm] ## adjacency matrix with permuted edge
}


## Function to change equal edges between two matrices
nomatch <- function(A, ## Original Matrix
                    Ap, ## Matrix to change
                    n.change){ ## Number of edges to change
  
  if(!isSymmetric(A) | !isSymmetric(Ap)) stop("Must have symmetric matrices")
  if(nrow(A) != nrow(Ap) | ncol(A) != ncol(Ap)) stop("Must have same size matrices")
  
  changed <- 0
  ## for loops are bad but this should be fast
  for(i in sample(1:nrow(A))){
    for(j in sample(1:nrow(A))){
      if(changed == n.change) break
      
      if(i != j){ ## Not dealing with diagonal
        if(A[i,j] == 0 & Ap[i,j] == 0){
          Ap[i,j] <- Ap[j,i] <- 1
          changed <- changed + 1
          break
        }
      }
    }
  }
  
  Ap[A == 1 & Ap == 1] <- 0
  Ap
}

## Function to make adj matrix more dense
increase.edge.density <- function(adj, edge.prob){
  for(j in 2:nrow(adj)){
    for(i in 1:(j-1)){ ## 1:(j-1) only deals with lower triangle
      if(adj[j,i] == 0){
        adj[j,i] <- adj[i,j] <- rbinom(1,1,edge.prob)  ## X% chance to add new edge
      } 
    }
  }
  adj
}


# Old stuff ---------------------------------------------------------------

# * 2 Moment Score --------------------------------------------------------

## Scaled Chisq test
ScaleChi <- function(P, K){
  ## Pieces of Itilde
  Itt <- (tr(P%*%K%*% P%*%K))/2
  Its <- tr(P%*%K%*%P)/2
  Iss <- tr(P%*%P)/2
  
  Itilde <- Itt - Its %*% solve(Iss) %*% t(Its)
  e <- tr(P%*%K)/2
  
  c(ka = Itilde/(2*e), ## Scale
    nu = (2*e^2)/Itilde ## Degrees of Freedom
  )
}

CC_Chisq_Score <- function(K, Y, mX, out.type = "C"){
  ## Projection matrix
  P <- diag(nrow(mX)) - mX%*%solve(t(mX)%*%mX)%*%t(mX) 
  R <- P%*%Y ## Residual
  
  s2 <- sum(R^2)/(nrow(mX) - ncol(mX)) ## MSE = SS/rdf
  sc <- ( t(R)%*%K%*%R )/(2*s2)
  s_chi <- ScaleChi(P, K)
  
  Q <- sc/s_chi['ka'] ## Scaled chisq stat
  pVal <- pchisq(Q, df = s_chi['nu'], lower.tail = FALSE)
  
  c(Q = Q, pVal = pVal, s_chi['ka'], s_chi['nu'])
}

perf_scor <- function(X, Y, Z, graph, delta){
  
  NL <- as.matrix( graph.laplacian(graph, normalized = T) )
  RNL <- solve(diag(p) + delta*NL)
  
  Z <- scale(Z); rho <- median(dist(Z))
  ZL <- Z %*% RNL
  
  K <- Gaussian_kernel(rho, ZL)
  CC_Chisq_Score(K, Y, X)
  
  # edge_density1 = edge_density(graph),
  # V.Y = V.Y, V.e = V.e, V.e.net = V.e.net,
  # total_edge_n = total_edge_n, pos_def = pd)
}

CCdavies.c <- function(H0.form, data = NULL, K){
  
  ## Null model
  m0 <- lm(H0.form, data = data)
  
  ## residuals, sigma^2, and model matrix for projection
  R <- m0$resid
  s2 <- summary(m0)$sigma^2
  mX <- model.matrix(H0.form, data = data)
  
  ## Projection matrix and observed stats
  P <- diag(nrow(mX)) - mX %*% (solve( t(mX)%*%mX ) %*% t(mX))
  PKP <- P %*% K %*% P
  q <- c(R %*% K %*% R)/s2
  
  ## Eigen values
  ee <- eigen(PKP - q * P, symmetric = T)$values
  ll <- ee[abs(ee) >= 1e-10]
  
  dd <- davies(q, lambda = ll)
  dd
}

CCdavies.d <- function(H0.form, data = NULL, K){
  mX <- model.matrix(H0.form, data)
  Y <- factor( eval(H0.form[[2]], data) )
  
  if (nlevels(Y) != 2) stop('The phenotype must be binary\n')
  
  m0 <- glm(H0.form, family = binomial, data = data)
  
  mu  <- m0$fitted.values
  res.wk <- m0$residuals
  res <- Y - mu
  w <- mu*(1-mu)
  sqrtw <- sqrt(w)
  
  adj <- sum((sqrtw * res.wk)^2) 
  DX12 <- sqrtw * mX
  
  qrX <- qr(DX12)
  Q <- qr.Q(qrX)
  Q <- Q[, 1:qrX$rank, drop=FALSE]
  
  P0 <- diag(nrow(X1)) - Q %*% t(Q)
  
  DKD <- tcrossprod(sqrtw) * K
  tQK <- t(Q) %*% DKD
  QtQK <- Q %*% tQK 
  PKP <- DKD - QtQK - t(QtQK) + Q %*% (tQK %*% Q) %*% t(Q)
  q <- as.numeric(res %*% K %*% res) / adj
  ee <- eigen(PKP, symmetric = T, only.values=T)  		
  lambda <- ee$values[abs(ee$values) >= tol]
}

compquad_pval <- function(AE, method){
  mms <- c("davies", "farebrother", "imhof", "liu")
  mm <- mms[grepl(method, mms)]
  
  if(mm == 'davies'){
    pvl <- CompQuadForm::davies(AE$Q, lambda = AE$ev)
  }
  
  if(mm == 'farebrother'){
    pvl <- CompQuadForm::farebrother(AE$Q, lambda = AE$ev)
  }
  
  if(mm == 'imhof'){
    pvl <- CompQuadForm::imhof(AE$Q, lambda = AE$ev)
  } 
  
  if(mm == "liu"){
    pvl <- CompQuadForm::liu(AE$Q, lambda = AE$ev)
  }
  
  pvl
}

## Simes Pvalue

pAdjSimes <- function(p){
  order(p)*p
}


## Test Functions ----
## Calculates Normalized Laplacian of metabolite pathway
getAdj <- function(target_compound, metabo2){
  cvnames <- metabo2[metabo2$KEGG %in% target_compound, c("KEGG", "BIOCHEMICAL") ]
  varnames <- cvnames$BIOCHEMICAL
  cnames <- cvnames[which(varnames %in% names(analysis)), "KEGG"]
  ncomp <- length(cnames)
  
  A <- matrix(-1, length(cnames), length(cnames))
  diag(A) <- 0
  for(i in 1:ncomp){
    r1 <- reactions[which(compoundReaction==cnames[i])]
    j <- 1
    while (j < i){
      r2 <- reactions[which(compoundReaction==cnames[j])]
      common <- intersect(r1, r2)
      
      if(length(common) > 0) {
        A[i, j] <- A[j, i] <- 1
      } else A[i, j] <- A[j, i] <- 0
      
      j <- j + 1
    }
  }
  
  if(sum(A)==0) return(-1)
  else{
    return(A)
  }
}

## Kernel test including network information through Reg Norm laplacian
kernel_RNL_test <- function(path_id, formula.H0, data,
                            metabs, tau = 1, kernel = "G",
                            adjusted = FALSE){
  
  target_compound <- names(keggGet(path_id)[[1]]$COMPOUND)
  metabo2 <- metabo[!is.na(metabo$KEGG) & !duplicated(metabo$KEGG), c("KEGG", "BIOCHEMICAL")]
  
  ## Returns normalized Laplacian
  A <- getAdj(target_compound, metabo2)
  if (all(A == -1)) return(-1)
  else{
    ## Graph info
    gg <- graph_from_adjacency_matrix(A, mode="undirected")
    NL <- laplacian_matrix(gg, normalized = T)
    
    ## metabolite measurements
    varnames <- metabo2[metabo2$KEGG %in% target_compound, "BIOCHEMICAL"]
    
    # print(varnames)
    # assign('vars', varnames, envir = globalenv())
    vv <- varnames[varnames %in% names(metabs)]
    ZZ <- scale(metabs[, vv])
    
    ## Number of disjoint pieces
    disj <- sum(round(eigen(NL)$values, 10) == 0)
    
    RNL <- solve(diag(nrow(NL)) + tau*NL)
    rho <- median(dist(ZZ))
    Z <- ZZ %*% RNL
    
    if(kernel == "G") K <- Gaussian_kernel(rho=rho, Z)
    if(kernel == "L") K <- plyKern(Z, pow = 1, rho = 0)
    
    sc <- SKAT.c(formula.H0, data = data, K=K, adjusted = adjusted)
    
    return(list(sc, disj = disj, tau = tau))
  }
}

kernel_NL_test <- function(path_id, formula.H0, data, 
                           tau = 1, metabs, kernel = "G"){
  
  target_compound <- names(keggGet(path_id)[[1]]$COMPOUND)
  metabo2 <- metabo[!is.na(metabo$KEGG) & !duplicated(metabo$KEGG), c("KEGG", "BIOCHEMICAL")]
  
  ## Returns normalized Laplacian
  A <- getAdj(target_compound, metabo2)
  if (all(A == -1)) return(-1)
  else{
    NL <- graph_from_adjacency_matrix(A, mode="undirected") %>%
      laplacian_matrix(normalized = T)
    
    varnames <- metabo2[metabo2$KEGG %in% target_compound, "BIOCHEMICAL"]
    ZZ <- scale(metabs[, varnames[varnames %in% names(metabs)]])
    
    ## Number of disjoint pieces
    disj <- sum(round(eigen(NL)$values, 10) == 0)
    
    rho <- median(dist(ZZ))
    Z <- as.matrix(ZZ %*% NL)
    
    if(kernel == "G") K <- Gaussian_kernel(rho, Z)
    if(kernel == "L") K <- plyKern(Z, pow = 1, rho = 0)
    
    sc <- SKAT.c(formula.H0, data = data, K=K)
    
    return(sc)
  }
}

## Kernel test without network information
No_Net_test <- function(path_id, formula.H0, data,
                        metabs, kernel = "G"){
  
  target_compound <- names(keggGet(path_id)[[1]]$COMPOUND)
  metabo2 <- metabo[!is.na(metabo$KEGG) & !duplicated(metabo$KEGG), c("KEGG", "BIOCHEMICAL")]
  varnames <- metabo2[metabo2$KEGG %in% target_compound, "BIOCHEMICAL"]
  
  ZZ <- scale(metabs[, varnames[varnames %in% names(metabs)]])
  
  ## So the output will match the test with network info
  A <- getAdj(target_compound, metabo2)
  if (all(A==-1)) return(NA)
  else{
    rho <- median(dist(ZZ)) 
    
    if(kernel == "G") K <- Gaussian_kernel(rho, ZZ)
    if(kernel == "L") K <- plyKern(ZZ, pow = 1, rho = 0)
    
    sc <- SKAT.c(formula.H0, data = data, K=K)
    return(sc)
  }
}

plot_graph <- function(path_id){
  target_compound <- names(keggGet(path_id)[[1]]$COMPOUND)
  metabo2 <- metabo[!is.na(metabo$KEGG) & !duplicated(metabo$KEGG), c("KEGG", "BIOCHEMICAL")]
  
  cvnames <- metabo2[metabo2$KEGG %in% target_compound, c("KEGG", "BIOCHEMICAL") ]
  varnames <- cvnames$BIOCHEMICAL
  cnames <- cvnames[which(varnames %in% names(analysis)), "KEGG"]
  ncomp <- length(cnames)
  
  A <- matrix(-1, length(cnames), length(cnames))
  diag(A) <- 0
  for(i in 1:ncomp){
    r1 <- reactions[which(compoundReaction==cnames[i])]
    j <- 1
    while (j < i){
      r2 <- reactions[which(compoundReaction==cnames[j])]
      common <- intersect(r1, r2)
      
      if(length(common) > 0) {
        A[i, j] <- A[j, i] <- 1
      } else A[i, j] <- A[j, i] <- 0
      
      j <- j + 1
    }
  }
  if(sum(A)==0) return(-1)
  gg <- graph_from_adjacency_matrix(A, mode="undirected")
  return(list(gg, path_id))
}

simes <- function(p){
  sp <- sort(p)
  ss <- length(sp) * sp / order(sp) 
  ss[order(order(p))]
}
