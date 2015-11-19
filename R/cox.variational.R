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
#' @param beta.start initial values of the fixed effect coefficients (used for iterative estimation scheme)
#' @param tol tolerance for judging convergence of the algorithm
#' @param verbose if \code{TRUE}, the algorithm prints verbose updates on its progess
#' @export
cox.variational <- function(y, X, S, wt, beta.start, tol=sqrt(.Machine$double.eps), verbose=TRUE) {
  r <- ncol(S)
  p <- ncol(X)
  n <- nrow(X)

  # Use a variational approximation with independence assumption as the initial estimate for iteration
  initial <- cox.variational.indep(y=y, X=X, S=S, wt=wt, beta.start=beta.start, tol=tol, verbose=verbose, hess=FALSE)

  # Interpret the initial estimates
  beta <- initial$beta
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
    result <- conjugate.gradient(logV=log(V)[indx], objective=likelihood.bound.logV, gradient=score.logV, y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M, tol=tol)
    # result <- conjugate.gradient(logV=result$par, objective=likelihood.bound.logV, gradient=score.logV.rev, y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M, tol=tol)

    # Recover the variance estimate of the random effects from the logged upper triangle
    logV <- result$par
    V <- matrix(0, r, r)
    V[indx] <- exp(logV)
    diagV <- diag(V)
    V <- V + t(V)
    diag(V) <- diagV

    # Holding V fixed, stimate M, beta, and tau by iteratively reweighted least squares
    if (verbose) cat("Holding V fixed to estimate M, tau, and beta")
    cholV <- as.matrix(t(chol(V)))
    v <- exp(VariationalVar(cholV, S) / 2)
    pseudoCovar <- rbind(cbind(X, S), cbind(Matrix(0, length(M), p),  sqrt(tau/2) * diag(r)))
    converged <- FALSE
    norm.old <- sum(beta^2)
    while(!converged) {
      z <- eta + (y / v - mu) / mu
      pseudodata <- c(z, rep(0, length(M)))
      pois.model <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt * mu * v, rep(1, r)))

      eta <- cbind(X, S) %*% pois.model$coefficients
      mu <- exp(eta)

      if (verbose) cat(".")
      norm.new <- sum(pois.model$coefficients^2)
      if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
      norm.old <- norm.new
    }
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
