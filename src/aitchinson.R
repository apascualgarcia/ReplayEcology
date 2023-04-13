###################################
## aitchinson.R
###################################
# This script aims to find an estimation of
# the correlation between two compositional variables
# using Aitchinson's transformation, following the
# strategy proposed in SparCC (Fierdman & Alm Plos Comp Biol 2012)
###################################
## Author: Dr. Alberto Pascual-García
## Date Created: 2022-12-19
## MIT Licensed.
## Contact: apascualgarcia.github.io
#################################### 
##
## INPUT: Two matrices A and B with samples in rows and in columns relative abundances
## OUTPUT: The cross-correlation between columns of A and columns of B
#####################################

aitchinson = function(A, B, pseudocount = 1e-07){
  #browser()
  Alog = log10(A + pseudocount)
  Blog = log10(B + pseudocount)
  N = dim(A)[2]
  M = dim(B)[2]
  D = N + M 
  # --- Compute an estimate of basis variances and total sum
  t_i = vector("numeric", length = N)
  t_j = vector("numeric", length = M)
  for(i in 1:N){
    t_i[i] = var(Alog[ ,i])
  }
  for(j in 1:M){
    t_j[j] = var(Blog[ ,j])
  }
  Tot = (sum(t_i) + sum(t_j))/(2*(D-1))
  # --- Estimate correlations
  Corr = matrix(0, nrow = N, ncol = M)
  for(i in 1:N){
    for(j in 1:M){
      diff = Alog[, i] - Blog[, j]
      t_ij = var(diff)
      w_i = sqrt((t_i[i] - Tot)/(D - 2))
      w_j = sqrt((t_j[j] - Tot)/(D - 2))
      Corr[i, j] = (w_i^2 + w_j^2 - t_ij)/(2*w_i*w_j)
    }
  }
  return(Corr)
}