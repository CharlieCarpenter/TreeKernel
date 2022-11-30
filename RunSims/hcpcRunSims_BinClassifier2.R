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

source("TreeKernel/BinClassifier/TwoVar/hcpcPower15_BinClassifier2.R")
         
res15 <- c(bonPow41=bonPow41, simPow41=simPow41,
           bonPow42=bonPow42, simPow42=simPow42,
           bonPow42G1=bonPow42G1, bonPow42G3=bonPow42G3,
           simPow42G1=simPow42G1, simPow42G3=simPow42G3,
           bonT140=bonT140, simT140=simT140)

pltDat <- cbind(res15, p=15)

source("TreeKernel/BinClassifier/TwoVar/hcpcPower30_BinClassifier2.R")

res30 <- c(bonPow41=bonPow41, simPow41=simPow41,
           bonPow42=bonPow42, simPow42=simPow42,
           bonPow42G1=bonPow42G1, bonPow42G3=bonPow42G3,
           simPow42G1=simPow42G1, simPow42G3=simPow42G3,
           bonT140=bonT140, simT140=simT140)

pltDat <- rbind(pltDat, cbind(res30, p=30))

source("TreeKernel/BinClassifier/TwoVar/hcpcPower45_BinClassifier2.R")

res45 <- c(bonPow41=bonPow41, simPow41=simPow41,
           bonPow42=bonPow42, simPow42=simPow42,
           bonPow42G1=bonPow42G1, bonPow42G3=bonPow42G3,
           simPow42G1=simPow42G1, simPow42G3=simPow42G3,
           bonT140=bonT140, simT140=simT140)

pltDat <- rbind(pltDat, cbind(res45, p=45))

pltDat %<>% as.data.frame %>% rownames_to_column(var="Sim") %>% 
  mutate(Group = case_when(grepl("G1",Sim) ~ "Group 1",
                           grepl("G3",Sim) ~ "Group 3"),
         Padj = case_when(grepl("bon",Sim) ~ "Bonferroni",
                          grepl("sim",Sim) ~ "Simes"),
         Clusters = case_when(grepl("4", Sim)~4,
                              TRUE ~ 2))

pltDat