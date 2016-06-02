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
#' @param beta.start starting values for iteration to estimate the fixed effect coefficients
#' @param tau.start initial value of the precision of the random effects
#' @param tol tolerance for judging convergence of the algorithm
#' @param verbose logical indicating whether to write detailed progress reports to standard output
#' @param hess logical indicating whether to compute the hessian after convergence (slow)
#' @param sd logical indicating whether to report the estimated standard errors of the coefficients
#'
#' @return list composed of these elements:
#'
#' \code{beta}: estimated vector of fixed effect regression coefficients
#'
#' \code{M}: estimated mean vector for the posterior of the spatial random effects at the converged value of the variational approximation
#'
#' \code{diagV}: vector of diagonal entries of the estimated covariance matrix for the posterior of the spatial random effects at convergence of the variational approximation
#'
#' \code{ltau}: estimated precision component for the spatial random effect
#'
#' \code{hessian}: estimated hessian matrix for \code{beta}, \code{M}, \code{log(diagV)} and \code{ltau} at convergence (estimated by call to \code{optim})
#'
#' \code{sd}: estimated standard errors of \code{beta}, \code{M}, and \code{ltau} at convergence (estimated by call to \code{sdreport})
#'
#' \code{object}: converged model object, returned by \code{nlminb}
#'
#' \code{neg.log.lik}: negative of the variational lower bound on the marginal log-likelihood at convergence
#'
cox.variational.indep <- function(y, X, S, wt, beta.start, tau.start=100, tol=sqrt(.Machine$double.eps), verbose, hess, sd) {
  # Start by estimating an optimal log(tau), assuming the given beta.start and u=rep(0,p)
  beta <- beta.start
  tau <- tau.start
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

  # Calculate an initial estimate for V
  diagV <- rep(1, r)
  log.diag.V <- log(diagV)

  # Estimate variance of the variational approximation
  obj <- MakeADFun(data=list(
    y=y,
    X=X,
    S=S,
    wt=wt
  ),
  parameters=list(
    logV = log.diag.V,
    M=u,
    beta=beta,
    ltau=ltau),
  DLL="tmb_diagonal")

  # run the variational approximation
  res <- nlminb(obj$par, obj$fn, obj$gr, obj$he, control=list(iter.max=10000, eval.max=20000))

  # extract parameter estimates from the result
  beta <- tail(res$par, p+1)[1:p]
  log.diag.V <- res$par[1:r]
  u <- res$par[r + 1:r]
  ltau <- tail(res$par, 1)

  # bundle up the return object
  out <- list(beta=beta, M=u, diagV=log.diag.V, ltau=ltau, neg.loglik=res$value, object=res)
  if (hess) out$hessian <- optimHess(res$par, fn=obj$fn, gr=obj$gr, y=y, X=X, S=S, wt=wt)
  if (sd) out$sd <- sdreport(obj)

  # return the results
  out
}

