###############################
##
## Project: MetaboGuru
##
## Purpose: Signal to Noise Ratio
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2022-05-24
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

# Helpful functions ----
`%nin%` <- Negate(`%in%`)
library(tidyverse); library(magrittr)

source("~/Documents/Research/Current/MetaboGuru/Carpenter/RCode/TreeKernels/HCPC/hcpcSims.R")

n <- 200; nsim <- 1000; p <- 30
tau <- 1; sd.y <- 1.3688 
b0 <- rep(0.1,5)

## Make 1,1 to include main effect
## Make 1,0 for only main effect 
A <- c(0,1)

# nGrp = 4,  b = 1,0,0,0 ------------------------------------------------
nGrp <- 4
set.seed(2345)
z1 <- rnorm(30,0,0.5)

set.seed(123)
dat.list <- lapply(1:nsim, function(x){
  grf <- sample_grg(p, 0.25)
  covs <- data.frame(X1 = c(rpois(n/4, 0),
                            rpois(n/4, 5),
                            rpois(n/4, 10),
                            rpois(n/4, 20)),
                     X2 = runif(n, -0.15, 0.15),
                     X3 = runif(n, -1, 1),
                     X4 = rnorm(n, sd=5))
  
  list(grf=grf, covs=covs)
})
betas <- lapply(seq(0.1,2, length.out = 10), ## Tuning parameter to multiply betas
                function(a) a*z1)

sig.noise.data <- lapply(betas, function(b){
  
  zz <- list(b, rep(0,p), rep(0,p), rep(0,p))
  
  ## Generating data set from every graph
  Y.Z <- lapply(dat.list, function(dl){
    
    hcDat <- hcTreeIntDat(dl$grf, dl$covs, 
                          zz=zz, groups=nGrp,
                          n=n, p=p,
                          b0=b0, sd.y=sd.y)
                 
    list(hcDat=hcDat, Ratio = hcDat$VY/hcDat$Ve)
  })
  ## Pulling Ratio from every data set 
  (sig.noise <- Y.Z %>% map_dbl("Ratio") %>% mean)
  
  list(Y.Z = Y.Z, Sig.Noise = sig.noise)
})
(sn <- unlist(map(sig.noise.data, "Sig.Noise")))

sndRes <- plyr::ldply(sig.noise.data, function(snd){
  
  tst <- sapply(snd$Y.Z, function(sndYZ){
    hcd <- sndYZ$hcDat
    hct <- hcTreeTest(hcd, formula(Y~X1+X2+X3+X4), nGrp)
    
    min(simes(hct$pvalHC))
  })
  
  data.frame(minSimes=tst, sigNoise=snd$Sig.Noise)
})


# save(sndRes, file = "Carpenter/RCode/TreeKernels/HCPC/sndRes.RData")
# load("Carpenter/RCode/TreeKernels/HCPC/SimRes/sndRes.RData")

sndRes %>% 
  mutate(sigNoise.f = factor(round(sigNoise,2))) %>% 
  group_by(sigNoise.f) %>% 
  dplyr::summarise(Power = sum(minSimes<0.05, na.rm=T)/n(),
                   .groups='drop') %>% 
  ggplot(aes(sigNoise.f, Power, group = 1)) +
  geom_point(color="dodgerblue2") +
  geom_line(color="dodgerblue2") +
  theme_bw() +
  labs(title="Power curve from hierarchical clustering partitioning",
       x=TeX("$\\frac{Var(\\beta Z)}{Var(\\epsilon)}$")) +
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0.2,1,0.2),
                     labels = as.character(seq(0.2,1,0.2)))



