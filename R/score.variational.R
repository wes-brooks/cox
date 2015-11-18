score.fin <- function(par, y, X, S, wt, V) {
  p <- ncol(X)
  r <- ncol(S)

  beta <- par[1:p]
  M <- par[(1:r) + p]
  ltau <- tail(par, 1)

  #grad <- as.vector(colSums(sweep(S^2, 1, wt * mu * v / 2, '*')))
  eta <- as.vector(X %*% beta + S %*% M)
  mu <- exp(eta)
  tau <- exp(ltau)

  cholV <- chol(V)
  cholV <- t(as.matrix(cholV))
  v <- exp(VariationalVar(cholV, S) / 2)

  # Gradient w.r.t. beta:
  grad.beta <- as.vector(t(X) %*% Diagonal(x=wt) %*% (y - mu*v))

  # Gradient w.r.t. M:
  grad.M <- as.vector(t(S) %*% Diagonal(x=wt) %*% (y - mu*v)) - tau*M

  # Gradient w.r.t. ltau:
  grad.ltau <- -tau/2 *sum(M^2) + ncol(S)/2

  -c(grad.beta, grad.M, grad.ltau)
}



#--------------------------
# Gradient of the lower likelihood bound:
score <- function(y, X, S, beta, wt, ltau, M, V) {
  eta <- as.vector(X %*% beta + S %*% M)
  r <- ncol(S)
  tau <- exp(ltau)
  mu <- exp(eta)

  cholV <- chol(V)
  cholV <- t(as.matrix(cholV))
  v <- exp(VariationalVar(cholV, S) / 2)

  grad <- VariationalScore(mu, wt, tau, v, as.matrix(V), S)

  as.vector(grad)
}



score.V <- function(V, y, X, S, beta, wt, ltau, M) {
  score(y, X, S, beta, wt, ltau, M, V)
}



score.logV <- function(logV, y, X, S, beta, wt, ltau, M) {
  V <- matrix(0, ncol(S), ncol(S))
  indx <- which(!lower.tri(V))
  V[indx] <- exp(logV)
  diagV <- diag(V)
  V <- V + t(V)
  diag(V) <- diagV

  eta <- as.vector(X %*% beta + S %*% M)
  r <- ncol(S)
  tau <- exp(ltau)
  mu <- exp(eta)

  cholV <- chol(V)
  cholV <- t(as.matrix(cholV))
  v <- exp(VariationalVar(cholV, S) / 2)

  grad <- VariationalScoreLogV(mu, wt, tau, v, as.matrix(V), S)
  grad[indx]
}

score.logV.rev <- function(logV, y, X, S, beta, wt, ltau, M) {
  -score.logV(logV, y, X, S, beta, wt, ltau, M)
}
