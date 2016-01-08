#' Gradient of the marginal log-likelihood for a Cox process
#'
#' Calculate the gradient of the marginal log-likelihood with respect to coefficients beta and log variance component ltau. Uses the Laplace approximation.
#'
#'
gradient.laplace <- function(beta, ltau, y, X, S, wt, tol, verbose) {
  tryCatch( {
    tau <- exp(ltau)
    p <- ncol(X)
    r <- ncol(S)
    eta <- as.vector(X %*% beta)
    mu <- exp(eta)

    if (verbose) cat("Maximizing joint likelihood for u:")
    pseudoCovar <- rbind(S, sqrt(tau/2) * diag(ncol(S)))
    converged <- FALSE
    norm.old <- Inf
    while(!converged) {
      z <- as.vector(eta - X%*%beta + (y - mu) / mu)
      pseudodata <- c(z, rep(0, r))
      pois <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt*mu, rep(1, r)))

      eta <- as.vector(X %*% beta + S %*% pois$coefficients)
      mu <- exp(eta)

      norm.new <- sum(pois$coefficients^2)
      if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
      norm.old <- norm.new
      if (verbose) cat(".")
    }
    u <- pois$coefficients
    if (verbose) cat("done!\n")

    eta <- as.vector(X %*% beta + S %*% u)
    mu <- exp(eta)

    # Calculate the Hessian for the u vector:
    H <- as.matrix(t(S) %*% Diagonal(x=wt * mu) %*% S)
    H <- (t(H) + H) / 2
    diag(H) <- diag(H) + tau

    # Compute the Cholesky decomposition of the Hessian
    L <- chol(H)
    L <- t(L)
    # deriv <- DerLogDetCholBeta(L, S, X, beta, u, wt)

#     # Rate of change of H w.r.t. beta
#     DH.Dbeta <- array(data=0, dim=c(r,r,p))
#     for (j in 1:p) {
#       # Rate of change of H w.r.t. beta[j]:
#       DH.Dbeta[,,j] <- t(S) %*% sweep(S, 1, as.vector(mu*wt*X[,j]), '*')
#     }
#
#     # Rate of change of log-det-cholesky factor w.r.t. beta.
#     D2 <- rep(0, p)
#     for (j in 1:p) {
#       # Rate of change of log-det-cholesky factor w.r.t. beta[j]
#       D2[j] <- sum(deriv * DH.Dbeta[,,j])
#     }

    D2 <- DerLogDetCholBeta(L, S, X, beta, u, wt)

    # Rate of change of log-det-cholesky factor w.r.t. tau
    # Dtau <- sum(diag(deriv))
    Dtau <- DerLogDetCholTau(L) * tau

    -c(as.vector(t(X) %*% Diagonal(x=wt) %*% (y - exp(eta)) - D2[1:p]), -tau*sum(u^2)/2 + r/2 - Dtau)
  }, error=function(e) return(NA)
  )
}
