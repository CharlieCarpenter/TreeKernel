###############################
##
## Project: MetaboGuru
##
## Purpose: Power simulations for PCA competing method
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
source("/Volumes/CIDA/BRANCHES/KechrisUO1/CharlieCarpenter/RCode/TreeKernels/CompeteingMethods/CompetingFunctions.R")

# Data Generation ---------------------------------------------------------
n <- 200; nsim <- 2000; p <- 30

tau <- 1; sd.y <- 1.3688 
b0 <- 0.1

set.seed(234)
graphL <- lapply(1:nsim, sample_grg, nodes=p, radius=0.25)
covL <- lapply(1:nsim, function(x){
  data.frame(X1 = interaction(rdunif(n, 2), rdunif(n, 2)),
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
                                   b0=b0, sd.y=sd.y,
                                   classifier='discrete'))

formula.H0 <- formula(Y~1)

pvalLHc41 <- lapply(hcDatL41, pcaTreeTest,
                    formula.H0 = formula.H0,
                    groups = nGrp, 
                    classifier = 'discrete')

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
                                   b0=b0, sd.y=sd.y,
                                   classifier='discrete'))

pvalLHc42 <- lapply(hcDatL42, pcaTreeTest,
                    formula.H0 = formula.H0,
                    groups = nGrp, 
                    classifier = 'discrete')

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
                                   b0=b0, sd.y=sd.y,
                                   classifier='discrete'))

pvalLHc40 <- lapply(hcDatL40, pcaTreeTest, 
                    formula.H0 = formula.H0,
                    groups = nGrp, 
                    classifier = 'discrete')

hcPval40 <- lapply(pvalLHc40, function(pl) pl$pvalHC)

bonHcPval40 <- lapply(hcPval40, p.adjust, method = 'bonferroni')
simHcPval40 <- lapply(hcPval40, simes)

bonT140 <- sum(sapply(bonHcPval40, min)<0.05)/nsim
simT140 <- sum(sapply(simHcPval40, min)<0.05)/nsim
