#' @useDynLib cox
NULL

library(Rcpp)
library(RcppEigen)
Rcpp::sourceCpp('src/LogDetDerChol.cpp')

