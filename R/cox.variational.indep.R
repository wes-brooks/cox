#' Estimate the regression parameters of a Cox process model using the variational approximation and assuming independence of the random effects.
#'
#' \code{cox.variational.indep} uses a variational approximation to estimate the parameters of a Cox process regression model with spatial random effects.
#' For this function, the variational approximation to the posterior distribution of the spatial random effects is a multivariate normal with
#' diagonal covariance matrix.
#'
#' @param y vector of response data
#' @param X design matrix for fized effects
#' @param S design matrix for the spatial random effects
#' @param wt vector of observation weights
#' @param beta.start starting values for iteration to estimate the fixed effect coefficients
#' @param tol tolerance for judging convergence of the algorithm
#' @param verbose logical indicating whether to write detailed progress reports to standard output
#' @param hess logical indicating whether to compute the hessian after convergence (slow)
#'
#' @return list composed of these elements:
#'
#' \code{beta}: estimated vector of fixed effect regression coefficients
#'
#' \code{M}: estimated mean vector for the posterior of the spatial random effects at the converged value of the variational approximation
#'
#' \code{diagV}: vector of diagonal entries of the estimated covariance matrix for the posterior of the spatial random effects at convergence of the variational approximation
#'
#' \code{ltau}: estimated precision component for the spatial random effect
#'
#' \code{hessian}: estimated hessian matrix for \code{beta}, \code{M}, \code{log(diagV)} and \code{ltau} at convergence (estimated by call to \code{optim})
#'
#' \code{neg.log.lik}: negative of the variational lower bound on the marginal log-likelihood at convergence
#'
#' @export
cox.variational.indep <- function(y, X, S, wt, beta.start, tau.start=100, tol=sqrt(.Machine$double.eps), verbose=TRUE, hess=TRUE) {
  # Start by estimating an optimal log(tau), assuming the given beta.start and u=rep(0,p)
  beta <- beta.start
  tau <- tau.start
  r <- ncol(S)
  p <- ncol(X)
  n <- nrow(X)

  # Prepare to begin iteration
  pseudoCovar <- rbind(S, sqrt(tau/2) * diag(r))
  eta <- as.vector(X %*% beta)
  mu <- exp(eta)
  norm.old <- Inf

  # Iterate to estimate u
  if (verbose) cat("Making an initial estimate of ltau")
  converged <- FALSE
  while(!converged) {
    z <- as.vector(eta - X %*% beta  + (y - mu) / mu)
    pseudodata <- c(z, rep(0, r))
    f1 <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt*mu, rep(1, r)))
    u <- f1$coefficients

    eta <- as.vector(X %*% beta + S %*% u)
    mu <- exp(eta)

    norm.new <- sum(u^2)
    if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
    norm.old <- norm.new

    if (verbose) cat(".")
  }
  tau <- r / norm.new
  ltau <- log(tau)
  if (verbose) cat(paste("done!\n ltau=", round(ltau, 3), "\n", sep=''))

  # Initial estimates of M, eta, and mu:
  # StS <- apply(S, 1, function(z) sum(z^2))
  M <- u
  eta <- as.vector(X %*% beta + S %*% M)
  mu <- exp(eta)

  # Calculate an initial estimate for V
#   conv <- FALSE
#   ltau.old <- -Inf
#   while (!conv) {
#     der2 <- -sum(exp(eta + exp(-ltau)*StS/2) * wt* (exp(-ltau)*StS/2 + (exp(-ltau)*StS/2)^2)) - sum(M^2)/2 * exp(ltau)
#     der <- sum(exp(eta + exp(-ltau)*StS/2) * wt * exp(-ltau)*StS/2) - exp(ltau)/2 * sum(M^2)
#     ltau <- ltau - der/der2
#
#     if (abs(ltau - ltau.old) < tol * (tol+abs(ltau.old))) {
#       conv <- TRUE
#     } else ltau.old <- ltau
#   }
#   tau <- exp(ltau)
  # diagV <- rep(1/tau, r)
  diagV <- rep(1, r)
  log.diag.V <- log(diagV)

  # Iterate and maximize the variational approximation to the marginal log-likelihood
  lik.old <- Inf
  conv.outer <- FALSE
  while (!conv.outer) {

    if (verbose) cat("Estimating variances for the variational approximation... ")
    converged <- FALSE
    ll.old <- Inf
    while(!converged) {
      # Estimate variance of the variational approximation
      res <- conjugate.gradient(M=M, logV=log.diag.V, objective=likelihood.bound.diagV, gradient=score.diagV, y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, tol=tol)
      log.diag.V <- res$logV
      diagV <- exp(log.diag.V)
      if (verbose) cat("done!\n")

      # Estimate mean vector for the variational approximation
      # First, calculate weights, offsets, and pseudodata for weighted least squares
      v <- exp(VariationalVarIndep(diagV, S) / 2)
      m <- as.vector(X %*% beta)
      pseudoCovar <- S
      pseudoCovar <- rbind(S, sqrt(tau/2) * diag(r))
      z <- eta - m + (y / v - mu) / mu
      pseudodata <- c(z, rep(0, length(M)))

      # Run weighted least squares and interpret the output
      pois <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt * mu * v, rep(1, r)))
      M <- pois$coefficients
      eta <- m + as.vector(S %*% M)
      mu <- exp(eta)

      # Estimate tau
      ltau <- log(ncol(S)) - log(sum(M^2) + sum(diagV))
      tau <- exp(ltau)

      # Check for convergence
      ll <- likelihood.bound.indep(M, log.diag.V, ltau, y, X, S, beta, wt)
      if (verbose) cat(paste("Checking convergence:\n Negative log-likelihood = ", round(ll, 3), "\n Convergence criterion = ", round(abs(ll - ll.old) / (tol * (tol + abs(ll.old))), 3), "\n\n"))
      if (abs(ll - ll.old) < tol * (tol + abs(ll.old)) | ll > ll.old) {
        converged <- TRUE
      } else ll.old <- ll
    }

    # Hold the variational approximation fixed to estimate the fixed effect coefficients
    # Calculate offset in preparation for IRLS
    pseudoCovar <- X
    m <- as.vector(S %*% M)
    converged <- FALSE

    # Iterate the IRLS algorithm
    while(!converged) {

      # Calculate pseudodata and do weighted least squares
      z <- eta - m + (y / v - mu) / mu
      pseudodata <- z
      pois <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt * mu * v))
      eta <- m + as.vector(X %*% pois$coefficients)
      mu <- exp(eta)

      # Check for convergence
      norm.new <- sum(pois$coefficients^2)
      if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
      norm.old <- norm.new
      if (verbose) cat(".")
    }
    beta <- pois$coefficients
    if (verbose) cat(" done!\n")

    # Check for convergence
    lik <- likelihood.bound.indep(M, log.diag.V, ltau, y, X, S, beta, wt)
    if (verbose) cat(paste(" ltau=", round(ltau, 3), "\n beta=", paste(round(beta, 3), collapse=', '), '\n\n', sep=''))
    if (verbose) cat(paste("Checking convergence:\n Negative log-likelihood = ", round(lik, 3), "\n Convergence criterion = ", round(abs(lik - lik.old) / (tol * (tol + abs(lik.old))), 3), "\n\n"))
    if (abs(lik - lik.old) < tol * (tol + abs(lik.old)) | lik > lik.old) {
      conv.outer <- TRUE
    } else lik.old <- lik
  }

  # Before returning, compute the Hessian if requested
  out <- list(beta=beta, M=M, diagV=diagV, ltau=ltau, neg.loglik=lik)
  if (hess) out$hessian <- optimHess(c(beta, M, log.diag.V, ltau), fn=likelihood.bound.fin.indep, gr=score.fin.indep, y=y, X=X, S=S, wt=wt)

  # Return the results
  out
}

