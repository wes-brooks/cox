gradient.joint <- function(par, y, X, S, u, wt) {
  p <- ncol(X)
  beta <- par[1:p]
  ltau <- tail(par, 1)
  -gradient.laplace(beta, ltau, y, X, S, u, wt)
}
