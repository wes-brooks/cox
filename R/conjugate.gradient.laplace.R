#' Maximize the marginal log likelihood with respect to regression coefficients
#' beta and spatial random effect precision component tau.
#'
#' @param objective function to be minimized (i.e. the marginal likelihood)
#' @param gradient function to calculate the gradient of the objective
#' @param y response variable for the regression function
#' @param X design matrix of covariates for the systematic portion of the model
#' @param S spatial random effect design matrix for the stochastic portion of the model
#' @param beta vector of initial regression coefficients for the fixed effects
#' @param u vector of initial values of the random effect loadings (coefficients)
#' @param wt vector giving the weight of each observation
#' @param ltau initial value for log of precision of the random effect loadings, \code{u}
#' @param verbose indicates whether to write detailed progress reports to standard output
#' @param tol proportional change in likelihood smaller than \code{tol} indicates convergence
#'
#' @return List consisting of \code{par}, the value of parameters that minimize
#' the objective function, and \code{value}, the value of the minimized objective
conjugate.gradient.laplace <- function(objective, gradient, y, X, S, beta, u, wt, ltau, verbose=TRUE, tol=sqrt(.Machine$double.eps)) {
  #Initial parameters:
  n <- nrow(S)
  r <- ncol(S)

  finished <- FALSE
  par <- c(beta, ltau)
  f.new <- objective(par, y=y, X=X, S=S, u=u, wt=wt, tol=tol, verbose=verbose)
  f.old <- Inf
  check <- Inf

  while(!finished) {

    #These iterations restart conjugacy:
    converged = FALSE
    iter <- 0
    while (!converged && iter<100) {
      iter <- iter+1

      #Prepare to iterate conjugate gradient descent:
      f.outer <- f.old
      t <- 1
      conv.inner <- FALSE
      i <- 0
      while(f.new < f.old && !converged && !conv.inner && i<p) {
        i <- i+1

        dir.new <- gradient(par, y=y, X=X, S=S, u=u, wt=wt, tol=tol, verbose=verbose)
        dir.new <- dir.new / sqrt(sum(dir.new^2))

        #First iteration, ignore conjugacy - thereafter, use it.
        #par.step is the vector of the new step (in parameter space)
        if (i==1) {
          par.step <- dir.new
        } else {
          conj <- (sum(dir.new^2) + sum(dir.new * dir.old)) / sum(dir.old^2)
          par.step <- dir.new + conj * s.old
        }

        #Find the optimal step size
        #Backtracking: stop when the loss function is majorized
        f.proposed <- objective(par + t*par.step, y=y, X=X, S=S, u=u, wt=wt, tol=tol, verbose=verbose)
        condition <- (f.proposed > f.new - sum((t*par.step)*dir.new) - 1/(2*t)*sum((t*par.step)^2))[1]
        condition <- f.proposed > f.new
        while(condition && t > .Machine$double.eps) {
          t = 0.5*t
          f.proposed <- objective(par + t*par.step, y=y, X=X, S=S, u=u, wt=wt, tol=tol, verbose=verbose)
          condition = (f.proposed > f.new - sum((t*par.step)*dir.new) - 1/(2*t)*sum((t*par.step)^2))[1]
          condition <- f.proposed > f.new

          #This is the final stopping rule: t gets so small that 1/(2*t) is Inf
          if (is.na(condition)) {
            converged = TRUE
            condition = FALSE
          }
        }

        #Find the optimal step
        par.proposed <- par + t*par.step

        #Make t a little bigger so next iteration has option to make larger step:
        t = t / 0.5 / 0.5

        #save for next iteration:
        dir.old <- dir.new
        s.old <- par.step

        #Only save the new parameters if they've decreased the loss function
        if (f.proposed < f.old)
          par <- par.proposed
        f.old <- f.new
        f.new <- f.proposed

        if ((f.old - f.new) < tol * (tol+abs(f.old))) conv.inner = TRUE
      }

      if (verbose) cat(paste("Iteration: ", iter, "; Objective: ", round(f.new, 3), "; Inner iterations: ", i, "\n", sep=""))
      if (abs(f.outer - f.new) < tol * (tol+abs(f.outer))) converged = TRUE
    }

    finished <- TRUE
  }

  list(par=par, value=f.new)
}
