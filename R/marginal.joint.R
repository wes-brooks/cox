# #--------------------------
# # Calculate the (negative) marginal log-likelihod
# marginal.joint <- function(par, y, X, S, u, wt) {
#   p <- ncol(X)
#   beta <- par[1:p]
#   ltau <- tail(par, 1)
#   marginal.laplace(beta, ltau, y, X, S, u, wt)
# }
#--------------------------
# Calculate the (negative) marginal log-likelihod
marginal.joint <- function(par, y, X, S, u, wt, tol, verbose) {
  p <- ncol(X)
  beta <- par[1:p]
  ltau <- tail(par, 1)
  marginal.laplace(beta, ltau, y, X, S, u, wt, tol, verbose)
}


marginal.joint.quick <- function(par, y, X, S, u, wt, verbose) {
  p <- ncol(X)
  beta <- par[1:p]
  ltau <- tail(par, 1)
  marginal.laplace.quick(beta, ltau, y, X, S, u, wt, verbose)
}


