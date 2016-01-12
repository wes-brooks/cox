#' Calculate the (negative) marginal log-likelihod
#'
#' @param beta vector of fixed-effect coefficients
#' @param ltau log-precision of the random effects
#' @param y vector of Poisson-distributed responses for DWPR
#' @param X matrix of fixed-effect covariates
#' @param S matrix of random-effect covariates
#' @param wt vector of observation weights
#' @param tol relative change threshold for convergence
#' @param verbose output detailed information about the progress of the algorithm?
#'
likelihood.laplace <- function(beta, ltau, y, X, S, wt, tol, verbose) {
  tryCatch( {
    tau <- exp(ltau)
    p <- ncol(X)
    r <- ncol(S)
    eta <- as.vector(X %*% beta)
    mu <- exp(eta)

    if (verbose) cat("Maximizing joint likelihood for u:")

    # Calculate offsets and pseudocovariates for penalized-IRLS
    m <- eta
    pseudoCovar <- rbind(S, sqrt(tau/2) * diag(ncol(S)))
    converged <- FALSE
    norm.old <- Inf

    # Iterate P-IRLS to convergence
    while(!converged) {

      # Generate pseudodata and then do weighted least squares
      z <- eta - m + (y - mu) / mu
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

    # Check for convergence
    nll <- -sum(wt * (y*eta - exp(eta))) + tau/2 * sum(u^2) - r/2*ltau + sum(log(diag(L)))
    if (verbose) cat(paste('beta=', paste(round(beta, 3), collapse=', '), "; ltau=", round(ltau, 3), "; neg-log-likelihood=", round(nll,3), '\n', sep=''))
    nll
    }, error=function(e) return(Inf)
  )
}


#' Calculate the (negative) marginal log-likelihod
#'
#' @param par vector of fixed-effect coefficients with log-precision parameter of the random effects appended
#' @param y vector of Poisson-distributed responses for DWPR
#' @param X matrix of fixed-effect covariates
#' @param S matrix of random-effect covariates
#' @param wt vector of observation weights
#' @param tol relative change threshold for convergence
#' @param verbose output detailed information about the progress of the algorithm?
#'
likelihood.joint <- function(par, y, X, S, wt, tol, verbose) {
  p <- ncol(X)
  beta <- par[1:p]
  ltau <- tail(par, 1)
  likelihood.laplace(beta, ltau, y, X, S, wt, tol, verbose)
}
