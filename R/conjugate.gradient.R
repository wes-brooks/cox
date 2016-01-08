#' Maximize the marginal log likelihood with respect to the covariance matrix of the spatial random effects
#'
#' @param objective function to be minimized (i.e. the marginal likelihood)
#' @param gradient function to calculate the gradient of the objective
#' @param y response variable for the regression function
#' @param X design matrix of covariates for the systematic portion of the model
#' @param S spatial random effect design matrix for the stochastic portion of the model
#' @param beta vector of initial regression coefficients for the fixed effects
#' @param wt vector giving the weight of each observation
#' @param ltau initial value for logarithm of prior distribution on the precision of the random effect loadings, \code{u} (precision component)
#' @param M vector of assumed posterior means of the spatial random effects
#' @param logV initial value of logarithm of the entries of the random effect posterior covariance matrix
#' @param verbose indicates whether to write detailed progress reports to standard output
#' @param tol proportional change in likelihood smaller than \code{tol} indicates convergence
#'
#' @return List consisting of \code{par}, the value of parameters that minimize
#' the objective function, and \code{value}, the value of the minimized objective
conjugate.gradient <- function(objective, gradient, y, X, S, beta, wt, ltau, M, logV, verbose=TRUE, tol=sqrt(.Machine$double.eps)) {
  #Initial parameters:
  n <- nrow(S)
  r <- ncol(S)

  finished <- FALSE
  f.new <- objective(y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M, logV=logV)
  check <- Inf

  while(!finished) {

    #These iterations restart conjugacy:
    converged = FALSE
    iter <- 0
    while (!converged && iter<100) {
      iter <- iter+1

      #Prepare to iterate conjugate gradient descent:
      f.outer <- objective(y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M, logV=logV)
      f.old <- Inf
      t <- 1
      conv.inner <- FALSE
      i <- 0

      if (verbose) cat("Iterating to estimate V")
      while(f.new < f.old && !converged && !conv.inner && i<(r*(r+1)/2 - 1)) {
        i <- i+1

        dir.new <- gradient(y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M, logV=logV)
        dir.new <- dir.new / sqrt(sum(dir.new^2))

        #First iteration, ignore conjugacy - thereafter, use it.
        #step is the vector of the new step (in parameter space)
        if (i==1) {
          step <- dir.new
        } else {
          conj <- (sum(dir.new^2) + sum(dir.new * dir.old)) / sum(dir.old^2)
          step <- dir.new + conj * s.old
        }

        #M.step <- step[1:r]
        #logV.step <- tail(step, length(step) - r)
        logV.step <- step

        #Find the optimal step size
        #Backtracking: stop when the loss function is majorized
        #f.proposed <- objective(M=M + t*M.step, logV=logV + t*logV.step, y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau)
        f.proposed <- objective(logV=logV + t*logV.step, y=y, X=X, S=S, beta=beta, M=M, wt=wt, ltau=ltau)
        #condition <- (f.proposed > f.new - sum((t*step)*dir.new) - 1/(2*t)*sum((t*step)^2))[1]
        condition <- (f.proposed > f.new)
        while(condition && t > .Machine$double.eps) {
          t = 0.5*t
          #f.proposed <- objective(M=M + t*M.step, logV=logV + t*logV.step, y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau)
          f.proposed <- objective(logV=logV + t*logV.step, y=y, X=X, S=S, beta=beta, M=M, wt=wt, ltau=ltau)
          #condition = (f.proposed > f.new - sum((t*step)*dir.new) - 1/(2*t)*sum((t*step)^2))[1]
          condition <- (f.proposed > f.new)

          #This is the final stopping rule: t gets so small that 1/(2*t) is Inf
          if (is.na(condition)) {
            converged = TRUE
            condition = FALSE
          }
        }

        #Find the optimal step
        #M.proposed <- M + t*M.step
        logV.proposed <- logV + t*logV.step

        #Make t a little bigger so next iteration has option to make larger step:
        t = t / 0.5 / 0.5

        #save for next iteration:
        dir.old <- dir.new
        s.old <- step

        #Only save the new parameters if they've decreased the loss function
        if (f.proposed < f.old) {
          #M <- M.proposed
          logV <- logV.proposed
        }
        f.old <- f.new
        f.new <- f.proposed

        if (verbose) cat(".")
        if ((f.old - f.new) < tol * (tol+abs(f.old))) conv.inner = TRUE
      }

      if (verbose) cat("done!\n")
      if (verbose) cat(paste("Iteration: ", iter, "; Objective: ", round(f.new, 3), "; Inner iterations: ", i, "\n", sep=""))
      if ((f.outer - f.new) < tol * (tol+abs(f.outer))) converged = TRUE
    }

    finished <- TRUE
  }

  #list(M=M, logV=logV, value=f.new)
  list(logV=logV, value=f.new)
}
