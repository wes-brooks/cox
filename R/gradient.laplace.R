#-------------------------
# Calculate the gradient of the marginal log-likelihood
gradient.laplace <- function(beta, ltau, y, X, S, u, wt, tol, verbose) {
  tryCatch( {
    tau <- exp(ltau)
    p <- ncol(X)
    eta <- X %*% beta + S %*% u
    mu <- exp(eta)
    r <- ncol(S)

    # Calculate the Hessian for the u vector:
    H <- as.matrix(t(S) %*% Diagonal(x=wt * mu) %*% S)
    H <- (t(H) + H) / 2
    diag(H) <- diag(H) + tau

    # Compute the Cholesky decomposition of the Hessian
    L <- chol(H)
    L <- t(L)
    deriv <- DerLogDetChol(L)

    # Rate of change of H w.r.t. beta
    DH.Dbeta <- array(data=0, dim=c(r,r,p))
    for (j in 1:p) {
      # Rate of change of H w.r.t. beta[j]:
      DH.Dbeta[,,j] <- t(S) %*% sweep(S, 1, as.vector(mu*wt*X[,j]), '*')
    }

    # Rate of change of log-det-cholesky factor w.r.t. beta.
    D2 <- rep(0, p)
    for (j in 1:p) {
      # Rate of change of log-det-cholesky factor w.r.t. beta[j]
      D2[j] <- sum(deriv * DH.Dbeta[,,j])
    }

    # Rate of change of log-det-cholesky factor w.r.t. tau
    Dtau <- sum(diag(deriv))


    c(as.vector(t(X) %*% Diagonal(x=wt) %*% (y - exp(eta)) - D2[1:p]), -tau*sum(u^2)/2 + r/2 - Dtau)
  }, error=function(e) return(NA)
  )
}
