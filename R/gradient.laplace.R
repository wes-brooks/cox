#-------------------------
# Calculate the gradient of the marginal log-likelihood
gradient.laplace <- function(beta, ltau, y, X, S, u, wt) {
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

    deriv <- LogDetDerChol2(L, S, X, mu, wt, tau)

    c(as.vector(t(X) %*% Diagonal(x=wt) %*% (y - exp(eta)) - deriv[1:p]), -tau*sum(u^2)/2 + r/2 - tail(deriv, 1))
  }, error=function(e) return(NA)
  )
}
