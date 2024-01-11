## 1. Initiation
# This file contains the packages used in the scripts 3-6 (objects, analysis, 
  #plotting, and supplement).

rm(list=ls())

#install necessary packages
install.packages("bipartite")
install.packages("vegan") 
install.packages("lme4") 
install.packages("tidyverse") 
install.packages("ggplot2") 
install.packages("viridis") 
install.packages("evolvability") 
install.packages('stringr') 
install.packages('ape') 
install.packages('picante')
install.packages("devtools")
devtools::install_github("paulponcet/oak")

#load libraries
library(bipartite) 
library(vegan) 
library(lme4) 
library(tidyverse) 
library(ggplot2) 
library(viridis) 
library(evolvability) 
library(stringr) 
library(ape) 
library(picante)
library(oak)
