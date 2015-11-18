marginal.laplace.beta <- function(par, ltau, y, X, S, u, wt, tol, verbose) {
  p <- ncol(X)
  beta <- par
  marginal.laplace(beta, ltau, y, X, S, u, wt, tol, verbose)
}
