#' Estimate the regression parameters of a Cox process model using the variational approximation with unconstrained covariance for the random effects.
#'
#' \code{cox.variational} uses a variational approximation to estimate the parameters of a Cox process regression model with spatial random effects.
#' For this function, the variational approximation to the posterior distribution of the spatial random effects is a multivariate normal with general
#' covariance matrix.
#'
#' @param y vector of response data
#' @param X design matrix for fized effects
#' @param S design matrix for the spatial random effects
#' @param wt vector of observation weights
#' @param beta.start starting values for iteration to estimate the fixed effect coefficients
#' @param tau.start initial value of the precision of the random effects
#' @param tol tolerance for judging convergence of the algorithm
#' @param verbose logical indicating whether to write detailed progress reports to standard output
#' @param hess logical indicating whether to estimate the Hessian matrix
#' @param sd logical indicating whether to report the estimated standard errors of the coefficients
#'
#' @return list of results containing the following elements:
#'
#' \code{beta}: estimated vector of fixed effect regression coefficients
#'
#' \code{M}: estimated mean vector for the posterior of the spatial random effects at the converged value of the variational approximation
#'
#' \code{V}: estimated covariance matrix for the posterior of the spatial random effects at the converged value of the variational approximation
#'
#' \code{ltau}: estimated precision component for the spatial random effect
#'
#' \code{hessian}: estimated hessian matrix for \code{beta}, \code{M}, and \code{ltau} at convergence (estimated by call to \code{optim})
#'
#' \code{sd}: estimated standard errors of \code{beta}, \code{M}, and \code{ltau} at convergence (estimated by call to \code{sdreport})
#'
#' \code{object}: converged model object, returned by \code{nlminb}
#'
#' \code{neg.log.lik}: negative of the variational lower bound on the marginal log-likelihood at convergence
#'
cox.variational <- function(y, X, S, wt, beta.start, tau.start=100, tol=sqrt(.Machine$double.eps), verbose, hess) {
  r <- ncol(S)
  p <- ncol(X)
  n <- nrow(X)

  # Use a variational approximation with independence assumption as the initial estimate for iteration
  initial <- cox.variational.indep(y=y, X=X, S=S, wt=wt, beta.start=beta.start, tau.start=tau.start, tol=tol, verbose=verbose, hess=FALSE, sd=FALSE)

  # Interpret the initial estimates
  beta <- initial$beta
  V <- diag(initial$diagV)
  ltau <- initial$ltau
  tau <- exp(ltau)
  M <- initial$M
  eta <- as.vector(X %*% beta + S %*% M)
  mu <- exp(eta)

  # Estimate variance of the variational approximation
  obj <- MakeADFun(data=list(
    y=y,
    X=X,
    S=S,
    wt=wt
  ),
  parameters=list(
    V=V,
    M=M,
    beta=beta,
    ltau=ltau),
  DLL="tmb_general")

  # run the variational approximation
  res <- nlminb(obj$par, obj$fn, obj$gr, control=list(iter.max=10000, eval.max=20000))

  # extract parameter estimates from the result
  beta <- tail(res$par, p+1)[1:p]
  V <- res$par[1:r]
  M <- res$par[r + 1:r]
  ltau <- tail(res$par, 1)

  # package up the results into the return object
  out <- list(beta=beta, M=M, V=V, ltau=ltau, neg.loglik=res$value, object=res)
  if (hess) out$hessian <- optimHess(c(beta, M, ltau), fn=likelihood.bound.fin, gr=score.fin, y=y, X=X, S=S, V=V, wt=wt)
  if (sd) out$sd <- sdreport(obj)

  # return the results
  out
}
