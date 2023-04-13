functionArgsList <- function (x, val)
{
  # this is a modification of the function appendList.R to merge two
  # lists, being x the default parameters for a given function and
  # val the parameters provided by the user, and returning default parameters
  # and those provided by the user, prevailing the latter in case there is
  # a match.
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  for (v in names(val)) {
    if (v %in% xnames){
      x[[v]] = val[[v]]
    }else{
      name.tmp=`v`
      x = c(x, setNames(as.list(val[[v]]),v))
    }
  }
  x
}


