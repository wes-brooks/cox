#' cox: A package for estimating Cox process regression models
#'
#' The \pkg{cox} package is used to estimate Cox process regression models. The
#' Cox process is a kind of mixed-effect model for spatial point processes. In
#' particular it includes a systematic or fixed effect regression portion and a
#' stochastic random effect portion. This package uses a fixed rank spatial
#' random effect (Cressie and Johannesson, 2008) for the stochastic portion and
#' estimates the regression parameters via a Poisson generalized linear model,
#' making use of numerical quadrature as in Renner et al. (2015). Finally, the
#' estimation is based on maximizing the model likelihood, marginal to the
#' random effects. Thus, the random effects are integrated out of the model via
#' a variational approximation.
#'
#'
#' @docType package
#' @name cox
NULL


#' Cox process parameter estimation
#'
#' Estimate the regression parameters of a Cox process model using the variational approximation.
#'
#' \code{cox} uses a variational approximation to estimate the parameters of a Cox process regression model with spatial random effects.
#' The variational approximation to the posterior distribution of the spatial random effects is a multivariate normal.
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
#' @param diagV logical: should the variational approximation to the covariance matrix of the random effects be diagonal?
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
#' \code{hessian}: estimated hessian matrix for \code{beta}, \code{M}, and \code{ltau} at convergence (estimated by call to \code{optimHess})
#'
#' \code{sd}: estimated standard errors of \code{beta}, \code{M}, and \code{ltau} at convergence (estimated by call to \code{sdreport})
#'
#' \code{object}: converged model object, returned by \code{nlminb}
#'
#' \code{neg.log.lik}: negative of the variational lower bound on the marginal log-likelihood at convergence
#'
#' @details Estimate regression coefficients of a spatial Cox process via the variational approximation.
#'
#' @export
cox <- function(y, X, S, wt, beta.start, tau.start=100, tol=sqrt(.Machine$double.eps), verbose=TRUE, hess=FALSE, diagV=FALSE, sd=TRUE) {

  # Dispatch the function call to the appropriate method: either independent random effects or not.
  if (diagV) res <- cox.variational.indep(y=y, X=X, S=S, wt=wt, beta.start=beta.start, tau.start=tau.start, tol=tol, verbose=verbose, hess=hess, sd=sd)
  else res <- cox.variational(y=y, X=X, S=S, wt=wt, beta.start=beta.start, tau.start=tau.start, tol=tol, verbose=verbose, hess=hess, sd=sd)

  res
}
