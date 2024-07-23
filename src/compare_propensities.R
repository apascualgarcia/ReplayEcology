###################################
#     compare_propensities.R  
##################################
# Author: Alberto Pascual-García
# Copyright (c)  Alberto Pascual-García,  2024
# Web:  apascualgarcia.github.io
# 
# Date: 2024-04-26
# Script Name: compare_propensities.R  
# Script Description: This script takes a df in which
# the different columns describe the propensity of
# each ASV (row) for each class (column). Additional
# columns provide quantiles of bootstrapped propensities
# and potentially further information.
# INPUT:
#     P_summary_all: A df with propensities and quantiles
#       for the different classes in columns and ASVs in rows.
#     pos_prop: A list of vectors, one vector for each class. Access
#        to each vector in the list is provided by the name of the class, 
#        and these names are used in the output.
#        Each vector x has three values: 
#              x[1] = column in the above df of the observed propensity of ASV i for being in class k 
#              x[2] = lower value of the confidence interval propensity of ASV i for being in class k
#              x[3] = upper confidence interval
# VALUE:
#     * Additional columns in P_summary_all describing the highest propensity
#       for the different classes of each ASV and if it is significant.
#     * A list with the ASVs having a significant propensity for each class.
#########################

compare_propensities = function(P_summary_all, pos_prop){
  classes = names(pos_prop)
  #classes = c("1","2","mix") # debug
  Nclasses = length(classes)
  #pos_prop = list() # debug
  #pos_prop[[classes[1]]] = c(1,2,3)
  #pos_prop[[classes[2]]] = c(5,6,7)
  #pos_prop[[classes[3]]] = c(9,10,11)
  #browser()
  ASV_sig_prop = list()
  
  prop_vec = unlist(lapply(pos_prop,FUN = function(x){x[1]}))
  ci_low_vec  = unlist(lapply(pos_prop,FUN = function(x){x[2]}))
  ci_high_vec  = unlist(lapply(pos_prop,FUN = function(x){x[3]}))
  # take for each set of propensities the columns of the top two propensities

  max_prop.cols = apply(P_summary_all[,prop_vec],   
                        MARGIN = 1,
                        FUN = function(x){prop_vec[sort.int(x, 
                                                   decreasing = T, 
                                                   index.return = T)$ix[c(1,2)]]})
  rownames(max_prop.cols)= c("max_1","max_2") # I take both just to double check everything is ok
  ci_low.cols = apply(P_summary_all[,prop_vec],   
                      MARGIN = 1,
                      FUN = function(x){ci_low_vec[sort.int(x, 
                                                          decreasing = T, 
                                                          index.return = T)$ix[1]]})
  ci_high.cols = apply(P_summary_all[,prop_vec],   
                    MARGIN = 1,
                    FUN = function(x){ci_high_vec[sort.int(x, 
                                                          decreasing = T, 
                                                          index.return = T)$ix[2]]})

  prop_sig = c()
  prop_sig_class = c()
  high_class = c()
  for(i in 1:dim(P_summary_all)[1]){
    prop_sig[i] = F
    prop_sig_class[i] = "none"
    high_class[i] = "none"
    max_val = P_summary_all[i, max_prop.cols[1,i]]
    if((length(max_val) == 0)){
      next
    }else if((max_val <= 0)| (max_val == abs(Inf)) # discard pathologies
       | (is.nan(max_val)) ){
      next
    }
    prop_diff = P_summary_all[i,ci_low.cols[i]]-
      P_summary_all[i,ci_high.cols[i]]
    high_class[i] = names(prop_vec)[which(prop_vec == max_prop.cols[1,i])]
    if(prop_diff > 0){
      prop_sig[i] = T
      prop_sig_class[i] = high_class[i]
    }
  }

  P_summary_all = cbind(P_summary_all,high_class,prop_sig,prop_sig_class)
  summary.table = table(prop_sig_class)
  ASV.list = list()
  for(class in names(pos_prop)){
    ASV.list[[class]] = rownames(P_summary_all)[which(prop_sig_class == class)]
  }
  # The following code is redundant and, so far, corrupts the output table
  # for(i in 1:Nclasses){
  #   label_sig =  paste0("sig_",classes[i])
  #   x1 = pos_prop[[classes[i]]][1]
  #   x2 = pos_prop[[classes[i]]][2]
  #   x3 = pos_prop[[classes[i]]][3]
  #   k = 0
  #   label_mode_j = c()
  #   for(j in 1:Nclasses){
  #     if(i == j){next}
  #     k = k+1
  #     y1 = pos_prop[[classes[j]]][1]
  #     y2 = pos_prop[[classes[j]]][2]
  #     y3 = pos_prop[[classes[j]]][3]
  #     label_sig_part = paste0("sig_",classes[i],"to",classes[j])
  #     label_mode = paste0("mode_",classes[i],"to",classes[j])
  #     label_mode_j[k] = label_mode
  #     P_summary_all[,label_sig_part] = 
  #       apply(P_summary_all, MARGIN = 1,
  #                FUN = function(x){unlist(sigclass(x1,x2,x3,
  #                                                  y1,y2,y3)[[1]])})
  #     P_summary_all[, label_mode] = 
  #       apply(P_summary_all, MARGIN = 1,
  #                FUN = function(x){unlist(sigclass(x1,x2,x3,
  #                                                  y1,y2,y3)[[2]])})
  #   }
  #   id_sig = which((P_summary_all[,label_mode_j[1]] == 1)& 
  #                             (P_summary_all[,label_mode_j[2]] == 1))
  #   P_summary_all[,label_sig] = F
  #   P_summary_all[id_sig, label_sig] = T
  #   ASV_sig_prop[[classes[i]]] = rownames(P_summary_all)[id_sig] 
  # }
  return(list("P_summary" = P_summary_all, "ASV_sig" = ASV.list))
} 

