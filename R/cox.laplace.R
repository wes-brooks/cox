#' Estimate the regression parameters of a Cox process model using the Laplace approximation.
#'
#' \code{cox.laplace} uses a Laplace approximation to integrate spatial random effects out of a Cox process regression model.
#' The Laplace approximation works by maximizing the joint likelihood of the data and random effects, then approximating the
#' random effects as multivariate normal via a second-order Taylor expansion.
#'
#' @param y vector of response data
#' @param X design matrix for fized effects
#' @param S design matrix for the spatial random effects
#' @param wt vector of observation weights
#' @param beta.start initial values of the fixed effect coefficients (used for iterative estimation scheme)
#' @param tol tolerance for judging convergence of the algorithm
#' @param verbose if \code{TRUE}, the algorithm prints verbose updates on its progess
#' @export
cox.laplace <- function(y, X, S, wt, beta.start, tau.start=100, tol=sqrt(.Machine$double.eps), verbose=TRUE) {

  # Start by estimating an optimal log(tau), assuming the given beta.start and u=rep(0,p)
  beta <- beta.start
  tau <- tau.start
  ltau <- log(tau)
  r <- ncol(S)
  p <- ncol(X)
  n <- nrow(X)
  par <- c(beta, ltau)

  # Estimate the regression coefficients
  #res <- optim(beta, fn=marginal.laplace, gr=gradient.laplace, y=dat$resp, X=X, S=S, u=u, ltau=ltau, wt=wt, L=L, method="BFGS")
  # res <- optim(c(beta, ltau), fn=marginal.laplace, gr=gradient.laplace, y=y, X=X, S=S, u=u, wt=wt, method="BFGS")
  # res <- optim(par=ltau, fn=marginal.tau, beta=beta, y=y, X=X, S=S, u=u, wt=wt, method="BFGS")
  # res <- optim(par=c(beta, ltau), fn=marginal.joint, y=y, X=X, S=S, wt=wt, verbose=verbose, method="BFGS")
  #res <- optim(ltau, fn=marginal.ltau, y=y, X=X, S=S, beta=beta, wt=wt, verbose=verbose, tol=tol)
  #res <- optim(par, fn=marginal.joint, gr=gradient.joint, y=y, X=X, S=S, wt=wt, verbose=verbose, tol=tol, method="BFGS")
  res <- optim(par, fn=marginal.joint, y=y, X=X, S=S, wt=wt, verbose=verbose, tol=tol, method="BFGS")
  beta <- res$par[1:p]
  tau <- tail(res$par, 1)

#
#   ltau <- res$par
#   tau <- exp(ltau)
#
#   # Compute eta and mu in order to begin the IRLS algorithm
#   eta <- X %*% beta
#   mu <- exp(eta)
#
#   # For fixed ltau, estimate beta and u by IRLS
#   if (verbose) cat("Holding tau fixed to estimate u.")
#   pseudoCovar <- rbind(cbind(X, S), cbind(matrix(0, ncol(S), ncol(X)), sqrt(tau/2) * diag(r)))
#   converged <- FALSE
#   norm.old <- Inf
#   while(!converged) {
#     z <- as.vector(eta + (y - mu) / mu)
#     pseudodata <- c(z, rep(0, r))
#     pois <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt*mu, rep(1, r)))
#
#     eta <- as.vector(cbind(X, S) %*% pois$coefficients)
#     mu <- exp(eta)
#
#     norm.new <- sum(pois$coefficients^2)
#     if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
#     norm.old <- norm.new
#     if (verbose) cat(".")
#   }
#   u <- tail(pois$coefficients, ncol(S))
#   beta <- pois$coefficients[1:ncol(X)]
#   if (verbose) cat ("done!\n")
#   if (verbose) cat(paste(" ltau=", round(ltau, 3), "\n initial beta=", paste(round(beta, 3), collapse=', '), '\n\n', sep=''))
#
#   if (verbose) cat(paste("Finalized estimation of u, estimating beta.\n", sep=''))
#   res <- optim(beta, fn=marginal.laplace.beta, y=y, ltau=ltau, X=X, S=S, u=u, wt=wt, verbose=verbose, tol=tol)
#   beta <- res$par

#   if (verbose) cat(paste("Finalized estimation of tau=", round(tau, 2), ", now estimating beta and its covariance matrix.\n", sep=''))
#   res.fin <- optimHess(c(beta, ltau), fn=marginal.laplace, y=y, X=X, S=S, u=u, wt=wt, verbose=verbose, tol=tol)

  # Estimate the Hessian of beta and tau, holding u, beta, and tau fixed at their maximizers
  if (verbose) cat("Finalized estimation of beta and tau, now estimate the hessian.\n")
  #res.fin <- optimHess(par=c(beta, ltau), fn=marginal.joint, gr=gradient.laplace, y=y, X=X, S=S, u=u, wt=wt, tol=tol, verbose=FALSE)
  #res.fin <- optimHess(par=c(beta, ltau), fn=marginal.joint, gr=gradient.joint, y=y, X=X, S=S, wt=wt, tol=tol, verbose=FALSE)
  res.fin=2

  list(beta=beta, ltau=ltau, hessian=res.fin, neg.loglik=res$value, u=u, res=res)
}
