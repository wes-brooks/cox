gradient.tau <- function(ltau, beta, y, X, S, u, wt) {
  -tail(gradient.laplace(beta, ltau, y, X, S, u, wt), 1)
}
