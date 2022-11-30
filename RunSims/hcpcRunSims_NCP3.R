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

source("/Volumes/CIDA/BRANCHES/KechrisUO1/CharlieCarpenter/RCode/TreeKernels/HCPC/NCP3/hcpcPower15_ncp3.R")

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

source("TreeKernel/NCP3/hcpcPower30_ncp3.R")

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

source("TreeKernel/NCP3/hcpcPower45_ncp3.R")

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

pltDat %<>% as.data.frame %>% rownames_to_column(var="Sim") %>% 
  mutate(Group = case_when(grepl("G1",Sim) ~ "Group 1",
                           grepl("G3",Sim) ~ "Group 3"),
         Padj = case_when(grepl("bon",Sim) ~ "Bonferroni",
                          grepl("sim",Sim) ~ "Simes"),
         Clusters = case_when(grepl("4", Sim)~4,
                              TRUE ~ 2))

pltDat
