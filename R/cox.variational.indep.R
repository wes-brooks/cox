#----------------------
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
#' @param beta.start initial values of the fixed effect coefficients (used for iterative estimation scheme)
#' @param tol tolerance for judging convergence of the algorithm
#' @param verbose if \code{TRUE}, the algorithm prints verbose updates on its progess
#' @export
cox.variational.indep <- function(y, X, S, wt, beta.start, tol=sqrt(.Machine$double.eps), verbose=TRUE, hess=TRUE) {
  # Start by estimating an optimal log(tau), assuming the given beta.start and u=rep(0,p)
  beta <- beta.start
  tau <- 100
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

    # Estimate V (update the variational approximation and maximize it)
    #diagV <- newton.indep(y, eta, S, wt, ltau, diagV)
    if (verbose) cat("Estimating variances for the variational approximation... ")
    # res <- optim(log.diag.V, fn=likelihood.bound.diagV, gr=score.diagV, y=y, M=M, X=X, S=S, beta=beta, wt=wt, ltau=ltau, control=list(reltol=.Machine$double.eps), method="BFGS")
    res <- conjugate.gradient(log.diag.V, objective=likelihood.bound.diagV, gradient=score.diagV, y=y, M=M, X=X, S=S, beta=beta, wt=wt, ltau=ltau, tol=tol)
    log.diag.V <- res$par
    diagV <- exp(log.diag.V)
    if (verbose) cat("done!\n")

    # Estimate M, beta, and tau for fixed V
    if (verbose) cat("Maximizing variational likelihood for beta and M: ")
    v <- exp(VariationalVarIndep(diagV, S) / 2)
    pseudoCovar <- rbind(cbind(X, S), cbind(Matrix(0, length(M), p),  sqrt(tau/2) * diag(r)))
    converged <- FALSE
    while(!converged) {
      z <- eta + (y / v - mu) / mu
      pseudodata <- c(z, rep(0, length(M)))
      pois <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt * mu * v, rep(1, r)))

      eta <- cbind(X, S) %*% pois$coefficients
      mu <- exp(eta)

      norm.new <- sum(pois$coefficients^2)
      if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
      norm.old <- norm.new
      if (verbose) cat(".")
    }
    beta <- pois$coefficients[1:ncol(X)]
    M <- tail(pois$coefficients, ncol(S))
    ltau <- log(ncol(S)) - log(sum(M^2) + sum(diagV))
    tau <- exp(ltau)
    if (verbose) cat(" done!\n")

    # Check for convergence
    lik <- likelihood.bound.indep(M, log.diag.V, ltau, y, X, S, beta, wt)

    if (verbose) cat(paste(" ltau=", round(ltau, 3), "\n beta=", paste(round(beta, 3), collapse=', '), '\n\n', sep=''))
    if (verbose) cat(paste("Checking convergence:\n Negative log-likelihood = ", round(lik, 3), "\n Convergence criterion = ", round(abs(lik - lik.old) / (tol * (tol + abs(lik.old))), 3), "\n\n"))
    if (abs(lik - lik.old) < tol * (tol + abs(lik.old))) {
      conv.outer <- TRUE
    } else lik.old <- lik
  }

  out <- list(beta=beta, M=M, diagV=diagV, ltau=ltau)

  #Compute the Hessian at the converged parameters:
  if (hess) out$hessian <- optimHess(c(beta, M, log.diag.V, ltau), fn=likelihood.bound.fin.indep, gr=score.fin.indep, y=y, X=X, S=S, wt=wt)

  out
}

