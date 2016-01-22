#' Variatonal lower bound on likelihood (eq 3.1 of Ormerod and Wand, 2012)
#'
likelihood.bound.indep <- function(M, log.diagV, ltau, y, X, S, beta, wt, StS) {
  r <- ncol(S)

  eta <- as.vector(X %*% beta + S %*% M)
  fixed <- as.vector(X %*% beta)
  mu <- exp(eta)
  tau <- exp(ltau)
  diagV <- exp(log.diagV)
  d <- VariationalVarIndep(diagV, S)
  v <- exp(d / 2)

  result <- sum(wt * (y * eta - mu * v)) # Expected conditional log-likelihood
  #result <- result + (r*ltau - tau*(sum(M^2) + sum(diagV))) / 2 # Add expectation of the prior on the random effects
  result <- result + (r*ltau - tau*(as.vector(t(M) %*% StS %*% M) + sum(d))) / 2 # Add expectation of the prior on the random effects
  result <- result + r/2*(1 + log(2*pi)) + sum(log(diagV)) # Add the entropy term
  -result
}


likelihood.bound.diagV <- function(logV, M, ltau, y, X, S, beta, wt, StS) {
  likelihood.bound.indep(M, logV, ltau, y, X, S, beta, wt, StS)
}


likelihood.bound.variational <- function(M, logV, ltau, y, X, S, beta, wt, StS) {
  likelihood.bound.indep(M, logV, ltau, y, X, S, beta, wt, StS)
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
