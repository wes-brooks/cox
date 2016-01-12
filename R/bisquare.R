#' Calculate weights via the bisquare kernel
#'
#' @param d vector of distances from the center of the bisquare kernel
#' @param bw bandwidth parameter for the bisquare kernel
#' @return The weight given to each of the distances \code{d} by a bisquare function with the given bandwidth \code{bw},
#' which is zero if \code{d>bw} and  \code{(1-(d/bw)^2)^2} otherwise
#'
bisquare <- function(d, bw)
{
  indx <- which(d < bw)
  resp <- rep(0, length(d))
  resp[indx] <- (1 - (d[indx]/bw)^2)^2

  return(resp)
}
