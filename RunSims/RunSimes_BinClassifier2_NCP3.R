###############################
##
## Project: MetaboGuru
##
## Purpose: Running TreeKernel simulations
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2022-06-02
##
## ---------------------------
## Notes:
## ---------------------------

# Helpful functions
`%nin%` <- Negate(`%in%`)
library(tidyverse); library(magrittr)

## Data Read In ----

source("CompeteingMethods/BinClassifier/TwoVar/Simes15_BinClasifier2_NCP3.R")

res15 <- c(bonPow41=bonPow41, simPow41=simPow41,
           bonPow42=bonPow42, simPow42=simPow42,
           bonPow42G1=bonPow42G1, bonPow42G3=bonPow42G3,
           simPow42G1=simPow42G1, simPow42G3=simPow42G3,
           bonT140=bonT140, simT140=simT140)

pltDat <- cbind(res15, p=15)

source("CompeteingMethods/BinClassifier/TwoVar/Simes30_BinClasifier2_NCP3.R")

res30 <- c(bonPow41=bonPow41, simPow41=simPow41,
           bonPow42=bonPow42, simPow42=simPow42,
           bonPow42G1=bonPow42G1, bonPow42G3=bonPow42G3,
           simPow42G1=simPow42G1, simPow42G3=simPow42G3,
           bonT140=bonT140, simT140=simT140)

pltDat <- rbind(pltDat, cbind(res30, p=30))

source("CompeteingMethods/BinClassifier/TwoVar/Simes45_BinClasifier2_NCP3.R")

res45 <- c(bonPow41=bonPow41, simPow41=simPow41,
           bonPow42=bonPow42, simPow42=simPow42,
           bonPow42G1=bonPow42G1, bonPow42G3=bonPow42G3,
           simPow42G1=simPow42G1, simPow42G3=simPow42G3,
           bonT140=bonT140, simT140=simT140)

pltDat <- rbind(pltDat, cbind(res45, p=45))

pltDat