#' Bagged CV bandwidth selector
#' 
#' @param x Vector. Sample.
#' @param r Positive integer. Size of the subsamples.
#' @param s Positive integer. Number of subsamples.
#' @param h0 Positive real number. Range over which to minimize, left bound.
#' @param h1 Positive real number. Range over which to minimize, right bound.
#' @param nb Positive integer. Number of bins to use in the \code{bw.ucv} function.
#' @param ncores Positive integer. Number of cores with which to parallelize the computations.
#' 
#' @details
#' Bagged cross-validation bandwidth for kernel density estimation.
#' 
#' @return Bagged CV bandwidth.
#' 
#' @examples
#' set.seed(1)
#' x <- rnorm(10^6)
#' bagcv(x, 5000, 100, 0.01, 1, 5000, 2)
#' 
#' @export
bagcv <-
function(x,r,s,h0,h1,nb=r,ncores=parallel::detectCores())
{
  n <- length(x)
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  paroutput <- foreach::foreach(i=1:s,.combine=c) %dopar%{
    subx <- sample(x,size=r,replace=FALSE)
    return(stats::bw.ucv(subx,nb=nb,h0,h1))
  }
  parallel::stopCluster(cl)
  hmean <- mean(paroutput)*(r/n)^0.2
  return(hmean)
}
