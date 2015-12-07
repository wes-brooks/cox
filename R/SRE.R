#' Set up the multiresolution spatial basis functions
#'
#' @param loc matrix of locations (both point observations and quadrature locations) that are to be included in the Cox process model. The columns must be named \code{x} and \code{y}.
#' @param n.1 number of basis functions that make up the first (coarsest) resolution
#' @param n.res number of resolutions to use
#' @param bw.scale multiplier used to set the bandwidths of the basis functions. Bandwidth for the \code{k}th resolution is set to \code{bw.scale * min.dist}, where \code{min.dist} is the minimum distance between basis function centers in the \code{k}th resolution (default is 1.5).
#' @param res.scale relative increase in the number of basis functions to use for each resolution (default is 3)
#' @param xlim bounds (in the x dimension) of the region to be spanned by the basis functions (optional)
#' @param ylim bounds (in the y dimension) of the region to be spanned by the basis functions (optional)
#'
#' @return \code{SRE}, the spatial random effect design matrix, with a row for each location and a column for each basis function. The \code{(i,j)}th entry of this matrix is the value of basis function \code{j} at location \code{i}.
#'
#' @export
SRE <- function(loc, n.1=16, n.res=3, bw.scale=1.5, res.scale=3, xlim=NULL, ylim=NULL) {
  if (!is.null(xlim)) XX <- range(loc$x) else XX <- xlim
  if (!is.null(ylim)) YY <- range(loc$y) else YY <- ylim
  n.grid <- n.1
  SRE <- matrix(NA, nrow(loc), 0)

  if (n.res == 0) return("All possible resolution layers have been eliminated.\n")

  # Set up three different spatial resolutions:
  for (r in 1:n.res) {
    margin <- 1 / 2 / sqrt(n.grid)

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

  # Remove random effect components that don't project onto any locations
  indx <- which(colSums(SRE)==0)
  if (length(indx)>0) SRE <- SRE[,-indx]

  SRE
}
