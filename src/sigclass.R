###################################
#     sigclass.R  
##################################
# Author: Alberto Pascual-García
# Copyright (c)  Alberto Pascual-García,  2024
# Web:  apascualgarcia.github.io
# 
# Date: 2024-04-22
# Script Name: sigboot.R  
# Script Description: It simply verifies that a
# value x is significantly higher or lower than 
# another value y, given their lower (p1 and q1)
# and upper (p2 and q2) 95% quantiles. It returns
# sig = T if it is significant and mode = -1 if
# it is significantly lower and mode = 1 if it
# is significantly higher. (sig = F, mode = 0 otherw.)

sigclass = function(x, p1, p2, y, q1, q2){
  #browser()
  if((x == 0)| (x == abs(Inf)) | (is.nan(x)) | 
     (y == abs(Inf)) | (is.nan(y)) |
     (p2 == Inf) | q2 == Inf){ # pathological cases
    sig = F
    sign = 0
  }else{
    if(x < y){
      if(p2 >=  q1){
        sig = F
        sign = 0
      }else{ # x is significantly lower
        sig = T
        sign = -1
      }
    }else{
      if(p1 <=  q2){
        sig = F
        sign = 0
      }else{ # significantly higher
        sig = T
        sign = 1
      }
    }
  }
  return(list("significance" = sig,
              "mode" = sign))
}