#' Calculate the (negative) marginal log-likelihod
#'
marginal.joint <- function(par, y, X, S, wt, tol, verbose) {
  p <- ncol(X)
  beta <- par[1:p]
  ltau <- tail(par, 1)
  marginal.laplace(beta, ltau, y, X, S, wt, tol, verbose)
}


marginal.joint.quick <- function(par, y, X, S, wt, verbose) {
  p <- ncol(X)
  beta <- par[1:p]
  ltau <- tail(par, 1)
  marginal.laplace.quick(beta, ltau, y, X, S, wt, verbose)
}


gradient.joint <- function(par, y, X, S, wt, tol, verbose) {
  p <- ncol(X)
  beta <- par[1:p]
  ltau <- tail(par, 1)
  gradient.laplace(beta, ltau, y, X, S, wt, tol, verbose)
}
