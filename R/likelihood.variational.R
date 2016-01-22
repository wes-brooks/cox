#' Variatonal lower bound on likelihood (eq 3.1 of Ormerod and Wand, 2012)
#'
#'
likelihood.bound <- function(y, X, S, beta, wt, ltau, M, V, StS) {
  r <- ncol(S)

  eta <- as.vector(X %*% beta + S %*% M)
  fixed <- as.vector(X %*% beta)
  mu <- exp(eta)
  tau <- exp(ltau)
  tryCatch( {
    cholV <- chol(V)
    cholV <- t(as.matrix(cholV))
    d <- VariationalVar(cholV, S)
    v <- exp(d / 2)

    # Compute the trace of the covariance matrix t(S) %*% V %*% S:
    tr <- sum(d)

    result <- sum(wt * (y * eta - mu * v)) #Expectation of conditional density
    result <- result + (r*ltau - tau*(as.vector(t(M) %*% StS %*% M) + tr)) / 2 #Expectation of prior on random effects
    result <- result + r/2*(1 + log(2*pi)) + sum(log(diag(cholV))) #Add the entropy term
    return(-result)
  }, error=function(e) return(Inf)
  )
}


likelihood.bound.V <- function(V, y, X, S, beta, wt, ltau, M, StS) {
  likelihood.bound(y, X, S, beta, wt, ltau, M, V, StS)
}


likelihood.bound.logV <- function(logV, y, X, S, beta, wt, ltau, M, StS) {
  V <- matrix(0, ncol(S), ncol(S))
  indx <- which(!lower.tri(V))
  V[indx] <- exp(logV)
  diagV <- diag(V)
  V <- V + t(V)
  diag(V) <- diagV

  likelihood.bound(y, X, S, beta, wt, ltau, M, V, StS)
}


likelihood.bound.fin <- function(par, y, X, S, wt, V, StS) {
  p <- ncol(X)
  r <- ncol(S)

  beta <- par[1:p]
  M <- par[(1:r) + p]
  ltau <- tail(par, 1)

  likelihood.bound(y, X, S, beta, wt, ltau, M, V, StS)
}





