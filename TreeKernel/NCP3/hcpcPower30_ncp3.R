###############################
##
## Project: MetaboGuru
##
## Purpose: Power simulations for HCPC
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

# Helpful functions
`%nin%` <- Negate(`%in%`)
library(tidyverse); library(magrittr)

source("/Volumes/CIDA/BRANCHES/KechrisUO1/CharlieCarpenter/RCode/TreeKernels/HCPC/hcpcSims.R")

# Data Generation ---------------------------------------------------------
n <- 200; nsim <- 2000; p <- 30

tau <- 1; sd.y <- 1.3688 
b0 <- rep(0.1,5)

set.seed(234)
graphL <- lapply(1:nsim, sample_grg, nodes=p, radius=0.25)
covL <- lapply(1:nsim, function(x){
  data.frame(X1 = c(rpois(n/4, 0),
                    rpois(n/4, 5),
                    rpois(n/4, 10),
                    rpois(n/4, 20)),
             X2 = runif(n, -0.15, 0.15),
             X3 = runif(n, -1,1),
             X4 = rnorm(n, sd=5))
})

# Group 4 -----------------------------------------------------------------

# * 1 active group --------------------------------------------------------

nGrp <- 4
set.seed(234)
zz <- list(rnorm(p,0,0.5), rep(0,p),
           rep(0,p), rep(0,p))
           
hcDatL41 <- mapply(hcTreeIntDat, graphL, covL,
                   SIMPLIFY = FALSE,
                   MoreArgs = list(zz=zz, groups=nGrp,
                                   n=n, p=p,
                                   b0=b0, sd.y=sd.y))

formula.H0 <- formula(Y~1)

pvalLHc41 <- lapply(hcDatL41, hcTreeTest, 
                    formula.H0 = formula.H0,
                    groups=nGrp, ncp=3)

hcPval41 <- lapply(pvalLHc41, function(pl) pl$pvalHC)

bonHcPval41 <- lapply(hcPval41, p.adjust, method = 'bonferroni')
simHcPval41 <- lapply(hcPval41, simes)

bonPow41 <- sum(sapply(bonHcPval41, min)<0.05)/nsim
simPow41 <- sum(sapply(simHcPval41, min)<0.05)/nsim

# * 2 active group --------------------------------------------------------

set.seed(234)
zz <- list(rnorm(p,0,0.3), rep(0,p),
           rnorm(p,0,0.3), rep(0,p))

hcDatL42 <- mapply(hcTreeIntDat, graphL, covL,
                   SIMPLIFY = FALSE,
                   MoreArgs = list(zz=zz, groups=nGrp,
                                   n=n, p=p,
                                   b0=b0, sd.y=sd.y))

pvalLHc42 <- lapply(hcDatL42, hcTreeTest, 
                    formula.H0 = formula.H0,
                    groups=nGrp, ncp=3)

hcPval42 <- lapply(pvalLHc42, function(pl) pl$pvalHC)

bonHcPval42 <- lapply(hcPval42, p.adjust, method = 'bonferroni')
simHcPval42 <- lapply(hcPval42, simes)

bonPow42 <- sum(sapply(bonHcPval42, min)<0.05)/nsim
simPow42 <- sum(sapply(simHcPval42, min)<0.05)/nsim

bonHcPval42G1 <- bonHcPval42 %>% map_dbl(~ .[[1]])
bonPow42G1 <- sum(bonHcPval42G1<0.05)/nsim

bonHcPval42G3 <- bonHcPval42 %>% map_dbl(~ .[[3]])
bonPow42G3 <- sum(bonHcPval42G3<0.05)/nsim

simHcPval42G1 <- simHcPval42 %>% map_dbl(~ .[[1]])
simPow42G1 <- sum(simHcPval42G1<0.05)/nsim

simHcPval42G3 <- simHcPval42 %>% map_dbl(~ .[[3]])
simPow42G3 <- sum(simHcPval42G3<0.05)/nsim

# * Type 1 ----------------------------------------------------------------

set.seed(234)
zz <- list(rep(0,p), rep(0,p),
           rep(0,p), rep(0,p))

hcDatL40 <- mapply(hcTreeIntDat, graphL, covL,
                   SIMPLIFY = FALSE,
                   MoreArgs = list(zz=zz, groups=nGrp,
                                   n=n, p=p,
                                   b0=b0, sd.y=sd.y))

pvalLHc40 <- lapply(hcDatL40, hcTreeTest, 
                    formula.H0 = formula.H0,
                    groups=nGrp, ncp=3)

hcPval40 <- lapply(pvalLHc40, function(pl) pl$pvalHC)

bonHcPval40 <- lapply(hcPval40, p.adjust, method = 'bonferroni')
simHcPval40 <- lapply(hcPval40, simes)

bonT140 <- sum(sapply(bonHcPval40, min)<0.05)/nsim
simT140 <- sum(sapply(simHcPval40, min)<0.05)/nsim

# Group 3 -----------------------------------------------------------------

covL <- lapply(1:nsim, function(x){
  data.frame(X1 = c(rpois(65, 0),
                    rpois(65, 10),
                    rpois(70, 20)),
             X2 = runif(n, -0.15, 0.15),
             X3 = runif(n, -1, 1),
             X4 = rnorm(n, sd=5))
})

# * 1 active group --------------------------------------------------------

nGrp <- 3
set.seed(234)
zz <- list(rnorm(p,0,0.5), rep(0,p),
           rep(0,p))

hcDatL31 <- mapply(hcTreeIntDat, graphL, covL,
                   SIMPLIFY = FALSE,
                   MoreArgs = list(zz=zz, groups=nGrp,
                                   n=n, p=p,
                                   b0=b0, sd.y=sd.y))



pvalLHc31 <- lapply(hcDatL31, hcTreeTest, 
                    formula.H0 = formula.H0,
                    groups=nGrp, ncp=3)

hcPval31 <- lapply(pvalLHc31, function(pl) pl$pvalHC)

bonHcPval31 <- lapply(hcPval31, p.adjust, method = 'bonferroni')
simHcPval31 <- lapply(hcPval31, simes)

bonPow31 <- sum(sapply(bonHcPval31, min)<0.05)/nsim
simPow31 <- sum(sapply(simHcPval31, min)<0.05)/nsim

# * 2 active group --------------------------------------------------------

set.seed(234)
zz <- list(rnorm(p,0,0.3), rnorm(p,0,0.3), rep(0,p))

hcDatL32 <- mapply(hcTreeIntDat, graphL, covL,
                   SIMPLIFY = FALSE,
                   MoreArgs = list(zz=zz, groups=nGrp,
                                   n=n, p=p,
                                   b0=b0, sd.y=sd.y))



pvalLHc32 <- lapply(hcDatL32, hcTreeTest, 
                    formula.H0 = formula.H0,
                    groups=nGrp, ncp=3)

hcPval32 <- lapply(pvalLHc32, function(pl) pl$pvalHC)

bonHcPval32 <- lapply(hcPval32, p.adjust, method = 'bonferroni')
simHcPval32 <- lapply(hcPval32, simes)

bonPow32 <- sum(sapply(bonHcPval32, min)<0.05)/nsim
simPow32 <- sum(sapply(simHcPval32, min)<0.05)/nsim

bonHcPval32G1 <- bonHcPval32 %>% map_dbl(~ .[[1]])
bonPow32G1 <- sum(bonHcPval32G1<0.05)/nsim

bonHcPval32G3 <- bonHcPval32 %>% map_dbl(~ .[[3]])
bonPow32G3 <- sum(bonHcPval32G3<0.05)/nsim

simHcPval32G1 <- simHcPval32 %>% map_dbl(~ .[[1]])
simPow32G1 <- sum(simHcPval32G1<0.05)/nsim

simHcPval32G3 <- simHcPval32 %>% map_dbl(~ .[[3]])
simPow32G3 <- sum(simHcPval32G3<0.05)/nsim

# * Type 1 ----------------------------------------------------------------

set.seed(234)
zz <- list(rep(0,p),rep(0,p), rep(0,p))

hcDatL30 <- mapply(hcTreeIntDat, graphL, covL,
                   SIMPLIFY = FALSE,
                   MoreArgs = list(zz=zz, groups=nGrp,
                                   n=n, p=p,
                                   b0=b0, sd.y=sd.y))



pvalLHc30 <- lapply(hcDatL30, hcTreeTest, 
                    formula.H0 = formula.H0,
                    groups=nGrp, ncp=3)

hcPval30 <- lapply(pvalLHc30, function(pl) pl$pvalHC)

bonHcPval30 <- lapply(hcPval30, p.adjust, method = 'bonferroni')
simHcPval30 <- lapply(hcPval30, simes)

bonT130 <- sum(sapply(bonHcPval30, min)<0.05)/nsim
simT130 <- sum(sapply(simHcPval30, min)<0.05)/nsim


# Group 2 -----------------------------------------------------------------

covL <- lapply(1:nsim, function(x){
  data.frame(X1 = c(rpois(n/2, 0),
                    rpois(n/2, 10)),
             X2 = runif(n, -0.15, 0.15),
             X3 = runif(n, -1, 1),
             X4 = rnorm(n, sd=5))
})

# * 1 active group --------------------------------------------------------

nGrp <- 2
set.seed(234)
zz <- list(rnorm(p,0,0.5), rep(0,p))

hcDatL21 <- mapply(hcTreeIntDat, graphL, covL,
                   SIMPLIFY = FALSE,
                   MoreArgs = list(zz=zz, groups=nGrp,
                                   n=n, p=p,
                                   b0=b0, sd.y=sd.y))

pvalLHc21 <- lapply(hcDatL21, hcTreeTest, 
                    formula.H0 = formula.H0,
                    groups=nGrp, ncp=3)

hcPval21 <- lapply(pvalLHc21, function(pl) pl$pvalHC)

bonHcPval21 <- lapply(hcPval21, p.adjust, method = 'bonferroni')
simHcPval21 <- lapply(hcPval21, simes)

bonPow21 <- sum(sapply(bonHcPval21, min)<0.05)/nsim
simPow21 <- sum(sapply(simHcPval21, min)<0.05)/nsim

# * Type 1 ----------------------------------------------------------------

set.seed(234)
zz <- list(rep(0,p), rep(0,p))

hcDatL20 <- mapply(hcTreeIntDat, graphL, covL,
                   SIMPLIFY = FALSE,
                   MoreArgs = list(zz=zz, groups=nGrp,
                                   n=n, p=p,
                                   b0=b0, sd.y=sd.y))

pvalLHc20 <- lapply(hcDatL20, hcTreeTest, 
                    formula.H0 = formula.H0,
                    groups=nGrp, ncp=3)

hcPval20 <- lapply(pvalLHc20, function(pl) pl$pvalHC)

bonHcPval20 <- lapply(hcPval20, p.adjust, method = 'bonferroni')
simHcPval20 <- lapply(hcPval20, simes)

sum(sapply(bonHcPval20, min)<0.05)/nsim
sum(sapply(simHcPval20, min)<0.05)/nsim

bonT120 <- sum(sapply(bonHcPval20, min)<0.05)/nsim
simT120 <- sum(sapply(simHcPval20, min)<0.05)/nsim





