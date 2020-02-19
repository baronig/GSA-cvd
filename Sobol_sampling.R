# Generate a Monte Carlo sample using Sobol' low-discrepancy quasi-random sequences
# James Keirstead
# 3 February 2012
# https://gist.github.com/jkeirstead/1730440
# 
# Random sampling with R's standard methods is inefficient for Monte Carlo analysis as
# the sampled values do not cover the parameter space evenly.  This Gist allows users
# to create parameter samples using Sobol' sequences to get around this problem.

# makeMCSample
# Makes a Monte Carlo sample using Sobol' sequences
#
# Parameters:
#  n = number of samples to draw
#  vals = list describing distributions, each consisting of:
#         - a variable name 'var', 
#         - the distribution name, 'dist', e.g. 'unif' for Uniform (see R's distribution names)
#         - the distribution parameters, 'params' (names vary by distribution)
#
# Returns:
#   a data frame with an index column "n" and each distribution in a named column
#
# modified by Gabriele Baroni
# 10.02.2020
# added the scrambling option
# eliminated the first colum with number of samples

makeMCSample <- function(n, vals,n.scramb,n.seed) {
  # Packages to generate quasi-random sequences
  # and rearrange the data
  require(randtoolbox)
  require(plyr)

  # Generate a Sobol' sequence
  sob <- sobol(n, length(vals), scrambling = n.scramb, seed = n.seed)

  # Fill a matrix with the values
  # inverted from uniform values to
  # distributions of choice
  samp <- matrix(rep(0,n*(length(vals))), nrow=n)
  for (i in 1:length(vals)) {
    l <- vals[[i]]
    dist <- l$dist
    params <- l$params
    fname <- paste("q",dist,sep="")
    samp[,i] <- do.call(fname,c(list(p=sob[,i]),params))
  }

  # Convert matrix to data frame and add labels
  samp <- as.data.frame(samp)
  names(samp) <- c(laply(vals, function(l) l$var))
  return(samp)
}

