#' Conjugate gradient optimization
#'
#' Maximize the marginal log-likelihood with respect to the log covariance matrix of the spatial random effects.
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
#'
#' @details Uses the method of conjugate gradient to maximize the marginal log-likelihood of the Cox process with respect to the log covariance matrix of the random effects. Uses backtracking to determine the step size, reducing the step size by half until stepping causes a decrease in the deviance. In the past I had used majorization-minimization to estimate an optimal step size, but the current approach seems to converge more quickly. The maximum number of steps before restarting conjugacy is the number of free parameters. In this case, that means the number of elements in the upper triangle (including the diagonal), which is r-choose-2 (where r is the dimension of the covariance matrix). Use only the upper diagonal because the covariance matrix because the covariance matrix is symmetric. Account for symmetry by doubling the gradient for all off-diagonal entries. Final convergence is when the relative decrease in the deviance is less than \code{tol}.
#'
conjugate.gradient <- function(objective, gradient, y, X, S, beta, wt, ltau, M, logV, verbose=TRUE, tol=sqrt(.Machine$double.eps)) {

  # Initial parameters:
  n <- nrow(S)
  r <- ncol(S)

  # Starting values for iteration
  finished <- FALSE
  f.new <- objective(y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M, logV=logV)
  check <- Inf

  # Iterate conjugate gradient until likelihood stops improving
  while(!finished) {

    # These iterations restart conjugacy:
    converged = FALSE
    iter <- 0
    while (!converged && iter<100) {
      iter <- iter+1

      # Prepare to iterate conjugate gradient descent:
      f.outer <- objective(y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M, logV=logV)
      f.old <- Inf
      t <- 1
      conv.inner <- FALSE
      i <- 0

      if (verbose) cat("Iterating to estimate V")
      while(f.new < f.old && !converged && !conv.inner && i<(r*(r+1)/2 - 1)) {
        i <- i+1

        # Compute the gradient of the likelihood function
        dir.new <- gradient(y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M, logV=logV)
        dir.new <- dir.new / sqrt(sum(dir.new^2))

        # First iteration ignores conjugacy.
        # step is the vector of the new step (in parameter space)
        if (i==1) {
          step <- dir.new
        } else {
          conj <- (sum(dir.new^2) + sum(dir.new * dir.old)) / sum(dir.old^2)
          step <- dir.new + conj * s.old
        }

        logV.step <- step

        # Find the optimal step size via backtracking
        f.proposed <- objective(logV=logV + t*logV.step, y=y, X=X, S=S, beta=beta, M=M, wt=wt, ltau=ltau)
        condition <- (f.proposed > f.new)
        while(condition && t > .Machine$double.eps) {
          t = 0.5*t
          f.proposed <- objective(logV=logV + t*logV.step, y=y, X=X, S=S, beta=beta, M=M, wt=wt, ltau=ltau)
          condition <- (f.proposed > f.new)

          #This is the final stopping rule: t gets so small that 1/(2*t) is Inf
          if (is.na(condition)) {
            converged = TRUE
            condition = FALSE
          }
        }

        # Only save the new parameters if they've decreased the loss function
        if (f.proposed < f.old)
          logV <- logV + t*logV.step

        # Make t a little bigger so next iteration has option to make larger step:
        t = t / 0.5 / 0.5

        # What's new is old in the next iteration:
        dir.old <- dir.new
        s.old <- step
        f.old <- f.new
        f.new <- f.proposed

        # Check for convergence
        if (verbose) cat(".")
        if ((f.old - f.new) < tol * (tol+abs(f.old))) conv.inner = TRUE
      }

      # Check for final convergence
      if (verbose) cat("done!\n")
      if (verbose) cat(paste("Iteration: ", iter, "; Objective: ", round(f.new, 3), "; Inner iterations: ", i, "\n", sep=""))
      if ((f.outer - f.new) < tol * (tol+abs(f.outer))) converged = TRUE
    }

    finished <- TRUE
  }

  # Return the result
  list(logV=logV, value=f.new)
}
