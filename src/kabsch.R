#===============================================================================
# PDBencode 
#' Kabsch Algorithm in R code
#'
#' Aligns two sets of points via rotation and translation after
#'   solving the superpositioning problem in an Eigenvalue Ansatz.
#'
#' Given two sets of points, with one specified as the reference set,
#' the other set will be rotated so that the RMSD between the two is minimized.
#' The format of the matrix is that there should be one row for each of
#' n observations, and the number of columns, d, specifies the dimensionality
#' of the points. The point sets must be of equal size and with the same
#' ordering, i.e. point one of the second matrix is mapped to point one of
#' the reference matrix, point two of the second matrix is mapped to point two
#' of the reference matrix, and so on.
#'
#' @param P n x d matrix of reference points.
#' @param Q n x d matrix of points to align to to \code{pm}
#' @return Matrix \code{qm} rotated and translated so that the ith point
#'	is aligned to the ith point of \code{pm} in the least-squares sense.
#' @references
#' \url{https://en.wikipedia.org/wiki/Kabsch_algorithm}
#'
#' (C) 2019 Jamie Macpherson and Jens Kleinjung
#===============================================================================

kabsch_R = function(Q, P) {
  ## cast array into matrix
  #dim(Q) = c(3, 4); # APG we already provide a matrix
  Q.m = Q; # t(Q); # APG each point is a community, in rows, so no transposition
  
  #dim(P) = c(3, 4);
  P.m = P; #t(P);
  
  ## center objects
  # center Q.m
  Qc = scale(Q.m, center = TRUE, scale = FALSE)
  
  # center P
  Pc = scale(P.m, center = TRUE, scale = FALSE)
  
  ## determine the cross-covariance matrix
  C = crossprod(Pc, Qc)
  
  ## compute a single value decomposition of the
  ## cross-covariance matrix
  C.svd = svd(C)
  
  ## use the sign of the determinant to ensure a right-hand coordinate system
  ## (this is only required if the fragments are stereogenic systems, which
  ## they aren't)
  #d = det(C.svd$u) * det(C.svd$v)
  
  ## determine the optimal rotation matrix
  U = tcrossprod(C.svd$u, C.svd$v)
  
  ## rotate matrix P into Q
  Prot = Pc %*% U
  
  ## determine the RMSD between rotated P and Q
  # squared differences
  Dsr.sq = (Prot - Qc)**2
  
  rmsd = sqrt(sum(apply(Dsr.sq, 1, sum)) / dim(Dsr.sq)[1]);
  
  ## spit out the optimal fit RMSD between Q and P and the transformation matrix
  return(list("U" = U, "rmsd" = rmsd, 
            "C" = C,  "Prot" = Prot))
}