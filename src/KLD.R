#############################
# Kullback-Leibler divergence
#############################

KLD <- function(x,y) sum(x *log(x/y))