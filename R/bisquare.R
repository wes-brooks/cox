#--------------------------
# The bisquare weight function
bisquare <- function(d, bw)
{
  indx <- which(d < bw)
  resp <- rep(0, length(d))
  resp[indx] <- (1 - (d[indx]/bw)^2)^2

  return(resp)
}
