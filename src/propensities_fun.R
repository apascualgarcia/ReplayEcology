###################################
#     propensities_fun.R  
##################################
# Author: Alberto Pascual-García
# Copyright (c)  Alberto Pascual-García,  2024
# Web:  apascualgarcia.github.io
# 
# Date: 2024-04-18
# Script Name: propensities_fun.R  
# Script Description: This function takes a list with 
# ASVs tables (the whole table and subtables, one for each
# class identified in the whole set), and another list
# with the probability of a community for being in those
# classes, and it returns the propensity for a community
# of being in a class given that an ASV was observed with
# a given frequency in that class.
# The function was created to iterate the computation through
# bootstrapped realizations of the ASV tables.
#
#

propensities_fun = function(ASV.table.class.bin.list, P_class.list,
                            class_names){
  #browser()
  freq_ASV.list = list()
  Prop.list = list()
  
  # --- Compute the frequency in the whole dataset
  ASV.table.bin = ASV.table.class.bin.list[["all"]]
  n_ASV_all = colSums(ASV.table.bin)
  freq_ASV_all = n_ASV_all/dim(ASV.table.bin)[1]
  freq_ASV.list[["all"]] = freq_ASV_all  
  for(class in class_names){
    ASV.table.class.bin = ASV.table.class.bin.list[[class]]
    # extract frequency of each ASV in each class
    n_ASV_class = colSums(ASV.table.class.bin)
    freq_ASV_class = n_ASV_class/n_ASV_all #dim(ASV.table.class.bin)[1]
    freq_ASV.list[[class]] = freq_ASV_class
  }
  
  # --- Compute propensity for each taxon
  for(class in class_names){
    freq_ASV_class = freq_ASV.list[[class]]
    P_class = P_class.list[[class]]
    P = log(freq_ASV_class)-log(freq_ASV_all)-log(P_class)
    P[which(P == -Inf)] = 0
    Prop.list[[class]] = P
  }
  
  return(list("Prop.list" =Prop.list,
              "freq_ASV.list" = freq_ASV.list))
}