indexG1.check = function (x, cl, d = NULL, centrotypes = "centroids") 
{
  #browser()
  if (sum(c("centroids", "medoids") == centrotypes) == 0) 
    stop("Wrong centrotypes argument")
  if ("medoids" == centrotypes && is.null(d)) 
    stop("For argument centrotypes = 'medoids' d cannot be null")
  if (!is.null(d)) {
    if (!is.matrix(d)) {
      d <- as.matrix(d)
    }
    row.names(d) <- row.names(x)
  }
  n <- length(cl)
  k <- max(cl)
  if (is.null(dim(x))) {
    dim(x) <- c(length(x), 1)
  }
  centers <- matrix(nrow = k, ncol = ncol(x))
  for (i in 1:k) {
    x.k = x[cl == i, ]
    if (centrotypes == "centroids") {
      if (ncol(x) == 1) {
        centers[i, ] <- mean(x.k)
      }
      else {
        if (is.vector(x.k)) {
          centers[i, ] <- x.k
        }
        else {
          centers[i, ] <- apply(x.k, 2, mean)
        }
      }
    }
    else {
      centers[i, ] <- .medoid(x[cl == i, ], d[cl == i, 
                                              cl == i])
    }
  }
  if (centrotypes == "centroids") {
    allmean <- apply(x, 2, mean)
  }
  else {
    allmean <- .medoid(x, d)
  }
  dmean <- sweep(x, 2, allmean, "-")
  allmeandist <- sum(dmean^2) # APG --> old definition, would not be needed
  withins <- rep(0, k)
  x <- (x - centers[cl, ])^2
  bgss=0 # APG
  for (i in 1:k) {
    withins[i] <- sum(x[cl == i, ])
    nels=length(which(cl == i))
    bgss = bgss + nels * sum((centers[i,]-allmean)^2) # APG --> new definition (see doc package clusterCrit)
  }
  wgss <- sum(withins)
  #bgss <- allmeandist - wgss # APG --> old definition
  (bgss/(k - 1))/(wgss/(n - k))
}