# Delete all objects from memory:
rm(list=ls())

packages <- c("readxl")
sapply(packages, require,character.only=TRUE) #Note, the "character.only" argument is necessary here

loc.wd <- "C:/Users/jwambaug/git/Comptox-Caco2/Honda_Caco2"

load(paste0(loc.wd, "/r_data/processed/all_gut_data_25MAR2019.RData")) 

musther <- readxl(paste0(loc.wd, "/r_data/lit_data/Musther-2014-DTXSID.xlsx")) 