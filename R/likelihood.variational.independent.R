#-------------------------
# Variatonal lower bound on likelihood (eq 3.1 of Ormerod and Wand, 2012)
likelihood.bound.indep <- function(M, log.diagV, ltau, y, X, S, beta, wt) {
  r <- ncol(S)

  eta <- as.vector(X %*% beta + S %*% M)
  mu <- exp(eta)
  tau <- exp(ltau)
  diagV <- exp(log.diagV)
  v <- exp(VariationalVarIndep(diagV, S) / 2)

  result <- sum(wt * (y * eta - mu * v))
  result <- result + (ncol(S)*ltau + sum(log.diagV) - tau*(sum(M^2) + sum(diagV))) / 2
  result <- result + r/2*(1 + log(2*pi)) + sum(log(diagV)) / 2
  -result
}


likelihood.bound.diagV <- function(log.diagV, M, ltau, y, X, S, beta, wt) {
  likelihood.bound.indep(M, log.diagV, ltau, y, X, S, beta, wt)
}


likelihood.bound.fin.indep <- function(par, y, X, S, wt) {
  p <- ncol(X)
  r <- ncol(S)

  beta <- par[1:p]
  M <- par[(1:r)+p]
  log.diagV <- par[(1:r)+(r+p)]
  ltau <- tail(par, 1)

  likelihood.bound.indep(M, log.diagV, ltau, y, X, S, beta, wt)
}
