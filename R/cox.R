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
