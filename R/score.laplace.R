#' Gradient of the marginal log-likelihood for a Cox process
#'
#' Calculate the gradient of the marginal log-likelihood with respect to coefficients beta and log variance component ltau. Uses the Laplace approximation.
#'
#' @param beta fixed-effect coefficient vector
#' @param ltau log precision of the random effect coefficients
#' @param y vector of Poisson-distributed responses for DWPR
#' @param X matrix of fixed-effect covariates
#' @param S matrix of random-effect covariates
#' @param wt vector of observation weights
#' @param tol relative change threshold for convergence
#' @param verbose output detailed information about the progress of the algorithm?
#'
score.laplace <- function(beta, ltau, y, X, S, wt, tol, verbose) {
  tryCatch( {

    # Calculate some perliminaries
    tau <- exp(ltau)
    p <- ncol(X)
    r <- ncol(S)
    eta <- as.vector(X %*% beta)
    mu <- exp(eta)

    if (verbose) cat("Maximizing joint likelihood for u:")

    # Calculate offsets and pseudocovariates for Penalized-IRLS
    m <- eta
    pseudoCovar <- rbind(S, sqrt(tau/2) * diag(ncol(S)))
    converged <- FALSE
    norm.old <- Inf

    # Iterate P-IRLS to convergence
    while(!converged) {

      # Get pseudodata and do weighted least squares
      z <- as.vector(eta - m + (y - mu) / mu)
      pseudodata <- c(z, rep(0, r))
      pois <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt*mu, rep(1, r)))

      # Interpret the output
      u <- pois$coefficients
      eta <- m + as.vector(S %*% u)
      mu <- exp(eta)

      # Check for convergence
      norm.new <- sum(pois$coefficients^2)
      if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
      norm.old <- norm.new
      if (verbose) cat(".")
    }
    if (verbose) cat("done!\n")

    # Calculate the Hessian for the u vector:
    H <- as.matrix(t(S) %*% Diagonal(x=wt * mu) %*% S)
    H <- (t(H) + H) / 2
    diag(H) <- diag(H) + tau

    # Compute the Cholesky decomposition of the Hessian
    L <- chol(H)
    L <- t(L)

    # Rate of change of log-det-cholesky factor w.r.t. parameters
    dbeta <- DerLogDetCholBeta(L, S, X, beta, u, wt)
    dtau <- DerLogDetCholTau(L) * tau

    # Return the result
    -c(as.vector(t(X) %*% Diagonal(x=wt) %*% (y - exp(eta)) - dbeta), -tau*sum(u^2)/2 + r/2 - dtau)
  }, error=function(e) return(NA)
  )
}


#' Gradient of the marginal log-likelihood for a Cox process
#'
#' Calculate the gradient of the marginal log-likelihood with respect to coefficients beta and log variance component ltau. Uses the Laplace approximation.
#'
#' @param par vector of fixed-effect coefficients, with the log-precision of the random effects appended
#' @param y vector of Poisson-distributed responses for DWPR
#' @param X matrix of fixed-effect covariates
#' @param S matrix of random-effect covariates
#' @param wt vector of observation weights
#' @param tol relative change threshold for convergence
#' @param verbose output detailed information about the progress of the algorithm?
#'
score.joint <- function(par, y, X, S, wt, tol, verbose) {
  p <- ncol(X)
  beta <- par[1:p]
  ltau <- tail(par, 1)
  score.laplace(beta, ltau, y, X, S, wt, tol, verbose)
}

