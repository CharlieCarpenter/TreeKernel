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

source("CompeteingMethods/BinClassifier/OneVar/Simes15_BinClassifier.R")

res15 <- c(bonPow41=bonPow41, simPow41=simPow41,
           bonPow42=bonPow42, simPow42=simPow42,
           bonPow42G1=bonPow42G1, bonPow42G3=bonPow42G3,
           simPow42G1=simPow42G1, simPow42G3=simPow42G3,
           bonT140=bonT140, simT140=simT140,
           bonPow31=bonPow31, simPow31=simPow31,
           bonPow32=bonPow32, simPow32=simPow32,
           bonPow32G1=bonPow32G1, bonPow32G3=bonPow32G3,
           simPow32G1=simPow32G1, simPow32G3=simPow32G3,
           bonT130=bonT130, simT130=simT130,
           bonPow21=bonPow21, simPow21=simPow21,
           bonT120=bonT120, simT120=simT120)

pltDat <- cbind(res15, p=15)

source("CompeteingMethods/BinClassifier/OneVar/Simes30_BinClassifier.R")

res30 <- c(bonPow41=bonPow41, simPow41=simPow41,
           bonPow42=bonPow42, simPow42=simPow42,
           bonPow42G1=bonPow42G1, bonPow42G3=bonPow42G3,
           simPow42G1=simPow42G1, simPow42G3=simPow42G3,
           bonT140=bonT140, simT140=simT140,
           bonPow31=bonPow31, simPow31=simPow31,
           bonPow32=bonPow32, simPow32=simPow32,
           bonPow32G1=bonPow32G1, bonPow32G3=bonPow32G3,
           simPow32G1=simPow32G1, simPow32G3=simPow32G3,
           bonT130=bonT130, simT130=simT130,
           bonPow21=bonPow21, simPow21=simPow21,
           bonT120=bonT120, simT120=simT120)

pltDat <- rbind(pltDat, cbind(res30, p=30))

source("CompeteingMethods/BinClassifier/OneVar/Simes45_BinClassifier.R")

res45 <- c(bonPow41=bonPow41, simPow41=simPow41,
           bonPow42=bonPow42, simPow42=simPow42,
           bonPow42G1=bonPow42G1, bonPow42G3=bonPow42G3,
           simPow42G1=simPow42G1, simPow42G3=simPow42G3,
           bonT140=bonT140, simT140=simT140,
           bonPow31=bonPow31, simPow31=simPow31,
           bonPow32=bonPow32, simPow32=simPow32,
           bonPow32G1=bonPow32G1, bonPow32G3=bonPow32G3,
           simPow32G1=simPow32G1, simPow32G3=simPow32G3,
           bonT130=bonT130, simT130=simT130,
           bonPow21=bonPow21, simPow21=simPow21,
           bonT120=bonT120, simT120=simT120)

pltDat <- rbind(pltDat, cbind(res45, p=45))

pltDat