#---------------------------
# Set up the spatial basis functions
#' @export
SRE <- function(loc, n.1=16, n.res=3, bw.scale=1.5, res.scale=3) {
  XX <- range(loc$x)
  YY <- range(loc$y)
  n.grid <- n.1
  SRE <- matrix(NA, nrow(loc), 0)

  if (n.res == 0) stop("All possible resolution layers have been eliminated.\n")

  # Set up three different spatial resolutions:
  for (r in 1:n.res) {
    margin <- runif(n=1, min=0, max=0.1)

    # Find the limits of this grid scale
    xlim <- c(margin, 1-margin) * diff(XX) + min(XX)
    ylim <- c(margin, 1-margin) * diff(YY) + min(YY)

    # Find the grid points at this scale
    len.xx <- round(sqrt(n.grid / (diff(YY) / diff(XX))))
    len.yy <- round(sqrt(n.grid / (diff(XX) / diff(YY))))
    xx <- seq(xlim[1], xlim[2], len=len.xx)
    yy <- seq(ylim[1], ylim[2], len=len.yy)
    closest <- min(xx[2] - xx[1], yy[2] - yy[1])
    points <- data.frame(x=rep(xx, times=len.yy), y=rep(yy, each=len.xx))

    # Evaluate the basis functions at the locations given in loc
    dist <- sqrt(outer(loc$x, points$x, '-')^2 + outer(loc$y, points$y, '-')^2)
    SRE <- cbind(SRE, apply(dist, 2, function(z) bisquare(z, bw.scale*closest)))

    #Get the number of points for the next scale of the grid
    n.grid <- n.grid * res.scale
  }

  # Remove any random effect components that are orthgonal to the points
  indx <- which(colSums(SRE)==0)
  if (length(indx)>0) SRE <- SRE[,-indx]

  # Final resolution is used for the grid to estimate K and tau
  len.xx <- round(sqrt(n.grid / (diff(YY) / diff(XX))))
  len.yy <- round(sqrt(n.grid / (diff(XX) / diff(YY))))
  xlim <- XX
  ylim <- YY
  xx <- seq(xlim[1], xlim[2], len=len.xx)
  yy <- seq(ylim[1], ylim[2], len=len.yy)
  closest <- min(xx[2] - xx[1], yy[2] - yy[1])
  points <- data.frame(x=rep(xx, times=len.yy), y=rep(yy, each=len.xx))
  dist <- sqrt(outer(loc$x, points$x, '-')^2 + outer(loc$y, points$y, '-')^2)
  bins <- apply(dist, 2, function(z) ifelse(z < bw.scale*closest, 1, 0))

  # Remove any empty bins
  indx <- which(colSums(bins)==0)
  if (length(indx)>0) bins <- bins[,-indx]

  return(list(SRE=SRE, bins=bins))
}
