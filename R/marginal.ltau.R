# #--------------------------
# # Calculate the (negative) marginal log-likelihod
# marginal.tau <- function(ltau, beta, y, X, S, wt) {
#   marginal.laplace(beta, ltau, y, X, S, wt)
# }


#--------------------------
# Calculate the (negative) marginal log-likelihod
marginal.ltau <- function(ltau, beta, y, X, S, wt, verbose, tol=sqrt(.Machine$double.eps)) {
  tryCatch( {
    tau <- exp(ltau)
    p <- ncol(X)
    r <- ncol(S)
    eta <- X %*% beta
    mu <- exp(eta)

    # Fix ltau and beta, maximize the joint likelihood through u
    if (verbose) cat("Maximizing joint likelihood for u: ")
    pseudoCovar <- rbind(cbind(X, S), cbind(matrix(0, ncol(S), ncol(X)), sqrt(tau/2) * diag(r)))
    converged <- FALSE
    norm.old <- Inf
    while(!converged) {
      z <- as.vector(eta + (y - mu) / mu)
      pseudodata <- c(z, rep(0, r))
      pois <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt*mu, rep(1, r)))

      eta <- as.vector(cbind(X, S) %*% pois$coefficients)
      mu <- exp(eta)

      norm.new <- sum(pois$coefficients^2)
      if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
      norm.old <- norm.new
      if (verbose) cat(".")
    }
    u <- tail(pois$coefficients, ncol(S))
    beta <- pois$coefficients[1:ncol(X)]
    if (verbose) cat(" done!\n")

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


