#--------------------------
# Calculate the (negative) marginal log-likelihod
marginal.laplace <- function(beta, ltau, y, X, S, u, wt, tol, verbose) {
  tryCatch( {
    tau <- exp(ltau)
    p <- ncol(X)
    r <- ncol(S)
    eta <- as.vector(X %*% beta + S %*% u)
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

    # Calculate the Hessian for the u vector:
    H <- as.matrix(t(S) %*% Diagonal(x=wt * mu) %*% S)
    H <- (t(H) + H) / 2
    diag(H) <- diag(H) + tau

    # Compute the Cholesky decomposition of the Hessian
    L <- chol(H)
    L <- t(L)

    r <- ncol(S)
    -sum(wt * (y*eta - exp(eta))) + tau/2 * sum(u^2) - r/2*log(tau) + sum(log(diag(L)))
  }, error=function(e) return(Inf)
  )
}




#--------------------------
# Calculate the (negative) marginal log-likelihod
marginal.laplace.quick <- function(beta, ltau, y, X, S, u, wt, verbose) {
  tryCatch( {
    tau <- exp(ltau)
    p <- ncol(X)
    eta <- as.vector(X %*% beta + S %*% u)
    mu <- exp(eta)

    # Calculate the Hessian for the u vector:
    H <- as.matrix(t(S) %*% Diagonal(x=wt * mu) %*% S)
    H <- (t(H) + H) / 2
    diag(H) <- diag(H) + tau

    # Compute the Cholesky decomposition of the Hessian
    L <- chol(H)
    L <- t(L)

    r <- ncol(S)
    -sum(wt * (y*eta - exp(eta))) + tau/2 * sum(u^2) - r/2*log(tau) + sum(log(diag(L)))
  }, error=function(e) return(Inf)
  )
}

