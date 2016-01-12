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
  grad.beta <- as.vector(t(X) %*% Matrix::Diagonal(x=wt) %*% (y - mu*v))

  # Gradient w.r.t. M:
  grad.M <- as.vector(t(S) %*% Matrix::Diagonal(x=wt) %*% (y - mu*v)) - tau*M

  # Gradient w.r.t. ltau:
  grad.ltau <- -tau/2 * (sum(M^2) + sum(diag(V))) + r/2

  -c(grad.beta, grad.M, grad.ltau)
}



#' Gradient of the variational likelihood bound w.r.t. log of covariance matrix
#'
#' @param logV logarithm of vectorized upper triangle of the covariance matrix of random-effect coefficients for the variational approximation
#' @param y vector of Poisson-distributed responses for DWPR
#' @param X matrix of fixed-effect covariates
#' @param S matrix of random-effect covariates
#' @param beta fixed-effect coefficient vector
#' @param wt vector of observation weights
#' @param ltau log precision of the random effect coefficients
#' @param M vector of means of random-effect coefficients for the variational approximation
#'
score.logV <- function(logV, y, X, S, beta, wt, ltau, M) {
  # Recover the covariance matrix of random-effect coefficients from its log upper triangle
  V <- matrix(0, ncol(S), ncol(S))
  indx <- which(!lower.tri(V))
  V[indx] <- exp(logV)
  diagV <- diag(V)
  V <- V + t(V)
  diag(V) <- diagV

  # Get linear predictor and means of the response
  eta <- as.vector(X %*% beta + S %*% M)
  mu <- exp(eta)
  r <- ncol(S)
  tau <- exp(ltau)

  # v is the part of the Gaussian MGF that depends on the covariance matrix
  cholV <- chol(V)
  cholV <- t(as.matrix(cholV))
  v <- exp(VariationalVar(cholV, S) / 2)

  # The covariance matrix is symmetric, and we use only the upper-triangular part.
  # So off-diagonal entries should count double to account for the entry across the diagonal.
  symmetrizer <- matrix(2, r, r)
  diag(symmetrizer) <- 1

  # Calculate the gradient and return it
  grad <- VariationalScoreLogV(mu, wt, tau, v, as.matrix(V), S)
  grad <- grad + DerLogDetChol(cholV) * symmetrizer * V;
  grad[indx]
}
