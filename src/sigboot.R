###################################
#     sigboot.R  
##################################
# Author: Alberto Pascual-García
# Copyright (c)  Alberto Pascual-García,  2024
# Web:  apascualgarcia.github.io
# 
# Date: 2024-04-22
# Script Name: sigboot.R  
# Script Description: It simply verifies that a
# value x is significantly higher or lower than zero,
# given its 0.05 quantiles q1 and q2

sigboot = function(x, q1, q2){
  #browser()
  if((x == 0)| (x == abs(Inf))
     | (is.nan(x))){
    sig = 0
    return(sig)}
  if(x > 0){
    if(q1 <= 0){
      sig = 0
    }else{
      if((x < q1)| (x > q2)){ # this should never happen
        sig = 0
      }else{
        sig = 1
      }
    }
  }else{
    if(q2 >= 0){
      sig = 0
    }else{
      if((x < q1)| (x > q2)){ # this should never happen
        sig = 0
      }else{
        sig = -1
      }
    }
  }
  return(sig)
}