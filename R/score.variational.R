#' Gradient of the variational negative log-likelihood with respect to beta, M, and ltau
#'
#' Gives the gradient of the variational negative log-likelihood with respect to coefficients \code{beta}, precision component \code{ltau}, and mean vector of the variational approximation, \code{M}. Within the \pkg{cox} package, this function is used to calculate gradients within a call to \code{optim} after the variational approximation has converged, as part of the process of calculating the Hessian.
#'
#' @param par vector composed of \code{beta}, \code{M}, and \code{ltau}
#' @param y vector of response values
#' @param X matrix of fixed effect covariates
#' @param S design matrix for the spatial random effects
#' @param wt vector of observation weights
#' @param V converged covariance matrix for the random effects, estimated by a variational approximation
#'
#' @return vector of gradients of the variational negative log-likelihood with respect to \code{par}
#'
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
  grad.ltau <- -tau/2 * (sum(M^2) + sum(diag(V))) + r/2

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


#
# score.V <- function(V, y, X, S, beta, wt, ltau, M) {
#   score(y, X, S, beta, wt, ltau, M, V)
# }



score.logV <- function(logV, y, X, S, beta, wt, ltau, M) {
  V <- matrix(0, ncol(S), ncol(S))
  indx <- which(!lower.tri(V))
  V[indx] <- exp(logV)
  diagV <- diag(V)
  V <- V + t(V)
  diag(V) <- diagV

  eta <- as.vector(X %*% beta + S %*% M)
  mu <- exp(eta)
  r <- ncol(S)
  tau <- exp(ltau)

  cholV <- chol(V)
  cholV <- t(as.matrix(cholV))
  v <- exp(VariationalVar(cholV, S) / 2)

  # The covariance matrix is symmetric, and we use only the upper-triangular part.
  # So off-diagonal entries should count double to account for the entry across the diagonal.
  symmetrizer <- matrix(2, r, r)
  diag(symmetrizer) <- 1

  grad <- VariationalScoreLogV(mu, wt, tau, v, as.matrix(V), S)
  grad <- grad + DerLogDetChol(cholV) * symmetrizer * V;
  grad[indx]
}


score.va <- function(M, logV, y, X, S, beta, wt, ltau) {
  #r <- ncol(S)
  #M <- par[1:r]
  #logV <- tail(par, length(par) - r)

  V <- matrix(0, ncol(S), ncol(S))
  indx <- which(!lower.tri(V))
  V[indx] <- exp(logV)
  diagV <- diag(V)
  V <- V + t(V)
  diag(V) <- diagV

  eta <- as.vector(X %*% beta + S %*% M)
  mu <- exp(eta)
  tau <- exp(ltau)

  cholV <- chol(V)
  cholV <- t(as.matrix(cholV))
  v <- exp(VariationalVar(cholV, S) / 2)


  # Gradient w.r.t. M:
  grad.M <- as.vector(t(S) %*% Diagonal(x=wt) %*% (y - mu*v)) - tau*M

  # The covariance matrix is symmetric, and we use only the upper-triangular part.
  # So off-diagonal entries should count double to account for the entry across the diagonal.
  symmetrizer <- matrix(2, r, r)
  diag(symmetrizer) <- 1

  grad <- VariationalScoreLogV(mu, wt, tau, v, as.matrix(V), S)
  grad <- grad + DerLogDetChol(cholV) * symmetrizer * V;
  c(grad.M, grad[indx])
}


score.logV.rev <- function(logV, y, X, S, beta, wt, ltau, M) {
  -score.logV(logV, y, X, S, beta, wt, ltau, M)
}
