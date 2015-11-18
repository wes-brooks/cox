#--------------------------
# Conjugate gradient to estimate the covariance matrix V of the spatial random effects
conjugate.gradient <- function(objective, gradient, y, X, S, beta, wt, ltau, M, logV, verbose=TRUE, tol=sqrt(.Machine$double.eps)) {
  #Initial parameters:
  n <- nrow(S)
  r <- ncol(S)

  finished <- FALSE
  f.new <- objective(logV, y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M)
  check <- Inf

  while(!finished) {

    #These iterations restart conjugacy:
    converged = FALSE
    iter <- 0
    while (!converged && iter<100) {
      iter <- iter+1

      #Prepare to iterate conjugate gradient descent:
      f.outer <- objective(logV, y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M)
      f.old <- Inf
      t <- 1
      conv.inner <- FALSE
      i <- 0

      if (verbose) cat("Iterating to estimate V")
      while(f.new < f.old && !converged && !conv.inner && i<(r*(r+1)/2 - 1)) {
        i <- i+1

        dir.new <- gradient(logV, y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M)
        dir.new <- dir.new / sqrt(sum(dir.new^2))

        #First iteration, ignore conjugacy - thereafter, use it.
        #step is the vector of the new step (in parameter space)
        if (i==1) {
          step <- dir.new
        } else {
          conj <- (sum(dir.new^2) + sum(dir.new * dir.old)) / sum(dir.old^2)
          step <- dir.new + conj * s.old
        }

        #Find the optimal step size
        #Backtracking: stop when the loss function is majorized
        f.proposed <- objective(logV + t*step, y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M)
        # condition <- (f.proposed > f.new - sum((t*step)*dir.new) - 1/(2*t)*sum((t*step)^2))[1]
        condition <- (f.proposed > f.new)
        while(condition && t > .Machine$double.eps) {
          t = 0.5*t
          f.proposed <- objective(logV + t*step, y=y, X=X, S=S, beta=beta, wt=wt, ltau=ltau, M=M)
          # condition = (f.proposed > f.new - sum((t*step)*dir.new) - 1/(2*t)*sum((t*step)^2))[1]
          condition <- (f.proposed > f.new)

          #This is the final stopping rule: t gets so small that 1/(2*t) is Inf
          if (is.na(condition)) {
            converged = TRUE
            condition = FALSE
          }
        }

        #Find the optimal step
        logV.proposed <- logV + t*step

        #Make t a little bigger so next iteration has option to make larger step:
        t = t / 0.5 / 0.5

        #save for next iteration:
        dir.old <- dir.new
        s.old <- step

        #Only save the new parameters if they've decreased the loss function
        if (f.proposed < f.old) {
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

  list(par=logV, value=f.new)
}
