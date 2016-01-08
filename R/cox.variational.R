#' Estimate the regression parameters of a Cox process model using the variational approximation with unconstrained covariance for the random effects.
#'
#' \code{cox.variational} uses a variational approximation to estimate the parameters of a Cox process regression model with spatial random effects.
#' For this function, the variational approximation to the posterior distribution of the spatial random effects is a multivariate normal with general
#' covariance matrix.
#'
#' @param y vector of response data
#' @param X design matrix for fized effects
#' @param S design matrix for the spatial random effects
#' @param wt vector of observation weights
#' @param beta.start starting values for iteration to estimate the fixed effect coefficients
#' @param tol tolerance for judging convergence of the algorithm
#' @param verbose logical indicating whether to write detailed progress reports to standard output
#'
#' @return list of results containing the following elements:
#'
#' \code{beta}: estimated vector of fixed effect regression coefficients
#'
#' \code{M}: estimated mean vector for the posterior of the spatial random effects at the converged value of the variational approximation
#'
#' \code{V}: estimated covariance matrix for the posterior of the spatial random effects at the converged value of the variational approximation
#'
#' \code{ltau}: estimated precision component for the spatial random effect
#'
#' \code{hessian}: estimated hessian matrix for \code{beta}, \code{M}, and \code{ltau} at convergence (estimated by call to \code{optim})
#'
#' \code{neg.log.lik}: negative of the variational lower bound on the marginal log-likelihood at convergence
#'
#' @export
cox.variational <- function(y, X, S, wt, beta.start, tau.start=100, tol=sqrt(.Machine$double.eps), verbose=TRUE) {
  r <- ncol(S)
  p <- ncol(X)
  n <- nrow(X)

  # Use a variational approximation with independence assumption as the initial estimate for iteration
  initial <- cox.variational.indep(y=y, X=X, S=S, wt=wt, beta.start=beta.start, tau.start=tau.start, tol=tol, verbose=verbose, hess=FALSE)

  # Interpret the initial estimates
  beta <- initial$beta
  #V <- matrix(min(diagV/100), r, r)
  #diag(V) <- initial$diagV
  V <- diag(initial$diagV)
  ltau <- initial$ltau
  tau <- exp(ltau)
  M <- initial$M
  eta <- as.vector(X %*% beta + S %*% M)
  mu <- exp(eta)

  # Iterate to maximize the variational approximation to the marginal log-likelihood
  lik.old <- Inf
  conv.outer <- FALSE
  while (!conv.outer) {

    # Estimate V by conjugate gradient descent
    indx <- which(!lower.tri(V))

    result <- conjugate.gradient(logV=log(V)[indx], objective=likelihood.bound.logV, gradient=score.logV, y=y, X=X, S=S, beta=beta, M=M, wt=wt, ltau=ltau, tol=tol)
    #result <- conjugate.gradient(M=M, logV=log(V)[indx], objective=likelihood.bound.va, gradient=score.va, y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, tol=tol)
    #result <- optim(c(M, log(V)[indx]), fn=likelihood.bound.va, gr=score.va, y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, method="BFGS")

    # Recover the variance estimate of the random effects from the logged upper triangle
    #M <- result$M
    logV <- result$logV

    #M <- result$par[1:r]
    #logV <- tail(result$par, length(result$par) - r)
    V <- matrix(0, r, r)
    V[indx] <- exp(logV)
    diagV <- diag(V)
    V <- V + t(V)
    diag(V) <- diagV

    # Holding V fixed, estimate M, beta, and tau by iteratively reweighted least squares
    if (verbose) cat("Holding V fixed to estimate M, tau, and beta")
    cholV <- as.matrix(t(chol(V)))
    v <- exp(VariationalVar(cholV, S) / 2)
    pseudoCovar <- rbind(cbind(X, S), cbind(Matrix(0, length(M), p),  sqrt(tau/2) * diag(r)))
    #m <- as.vector(S %*% M)
    #pseudoCovar <- X
    converged <- FALSE
    norm.old <- sum(beta^2)
    while(!converged) {
      z <- eta + (y / v - mu) / mu
      pseudodata <- c(z, rep(0, length(M)))
      pois.model <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt * mu * v, rep(1, r)))
      eta <- cbind(X, S) %*% pois.model$coefficients

      #z <- eta - m + (y / v - mu) / mu
      #pseudodata <- z
      #pois <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt * mu * v))
      #eta <- m + as.vector(X %*% pois$coefficients)

      mu <- exp(eta)

      if (verbose) cat(".")
      norm.new <- sum(pois$coefficients^2)
      if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
      norm.old <- norm.new
    }
    #beta <- pois$coefficients
    beta <- pois.model$coefficients[1:p]
    M <- tail(pois.model$coefficients, r)
    ltau <- log(r) - log(sum(M^2) + sum(diag(V)))
    tau <- exp(ltau)
    if (verbose) cat(paste("done!\n beta = ", paste(round(beta, 3), collapse=", "), "\n ltau = ", round(ltau, 3), "\n", sep=''))

    # Check for convergence
    lik <- likelihood.bound(y, X, S, beta, wt, ltau, M, V)
    if (verbose)cat(paste("Checking final convergence criterion:\n likelihood=", round(lik, 3), "\n convergence criterion=", round(abs(lik - lik.old) / (tol * (tol + abs(lik.old))), 3), "\n"))
    if (abs(lik - lik.old) < tol * (tol + abs(lik.old))) {
      conv.outer <- TRUE
    } else lik.old <- lik
  }

  res.fin <- optimHess(c(beta, M, ltau), fn=likelihood.bound.fin, gr=score.fin, y=y, X=X, S=S, V=V, wt=wt)
  list(beta=beta, M=M, V=V, ltau=ltau, hessian=res.fin, neg.loglik=lik)
}
