#' Simulate a Gaussian random field a set of points
#'
#' @param intercept intercept for the regression portion of the model
#' @param beta coefficients of the regression portion of the model
#' @param rando.corr
#' @param rho correlation between the
#' @param n.q number of quadrature points
#' @param n
#' @param output.loc save output to this folder
#' @param reps number of times to repeat this simulation setting
#' @export
#'
sim.cox <- function(intercept, beta, rando.corr, rho, n.q, n, output.loc, reps) {
  # #------------------
  # # Parameters
  # intercept <- c(1,2,5,10)
  # beta <- c(1, -1)
  # rando.corr <- c(0, 0)
  # rho <- matrix(0, p, p)
  # rho[1,] <- rho[,1] <- 0.2
  # diag(rho) <- 1


  #--------------------
  # Generate the observation points
  inc <- rexp(n.q + 1)
  x <- cumsum(inc[1:n.q]) / sum(inc)
  y <- runif(n.q)
  loc <- data.frame(x=x, y=y)

  # Generate the quadrature points
  n <- round(sqrt(n.q))
  x <- rep(seq(0, 1, length=n), each=n)
  y <- rep(seq(0, 1, length=n), times=n)
  loc <- rbind(loc, data.frame(x=x, y=y))


  #--------------------
  # Simulate the random effect
  kr <- RFoptions()$krige
  kr$locmaxn <- n.q+1
  re <- GRF(loc=data.frame(x=loc$x, y=loc$y), tau=0.1, sigma=0.5, spConform=FALSE)




  #---------------------
  # Simulate the covariates
  p = 3
  X = matrix(0, nrow(loc), 0)
  fields <- list()
  for(j in 1:p) {
    d <- GRF(loc=data.frame(x=loc$x, y=loc$y), tau=0.25, sigma=1)
    fields[[j]] <- d
    X = cbind(X, d@data[[1]])
  }

  # X = as.matrix((x + y) / 4)

  #
  #
  # #--------------------
  # # Apply correlation between covariates, random effect:
  # rho <- matrix(0.4, 2, 2)
  # diag(rho) <- 1
  # Z <- cbind(X, re@data[[1]]) %*% chol(rho)
  #




  obs.indx <- 1:n.q
  quad.indx <- (n.q+1):nrow(loc)
  # eta <- 2 + Z[obs.indx,1] + Z[obs.indx,2] #+ X[,1]
  eta <- as.vector(intercept + X %*% beta  + re@data[[1]])
  sup <- exp(max(eta))


  #-----------------------
  # Apply the simulated regression model and then resample the points to get the Cox process:

  count <- 0
  while(TRUE) {
    count <- count + 1

    # Grab some points for the homogeneous Poisson process
    n.ppp <- rpois(1, sup)
    indx.ppp <- sample(obs.indx, n.ppp)
    X.ppp <- X[indx.ppp,]

    write.csv(loc[indx.ppp,], file=paste("~/Dropbox/confounding/output/loc.ppp", count, "out", sep="."), row.names=FALSE)

    # Thin the points we just selected
    prob <- exp(eta[indx.ppp]) / sup
    cox.indx <- which(rbinom(n.ppp, size=1, prob=prob)==1)
    X.cox <- as.matrix(X.ppp[cox.indx,])
    n.cox <- length(cox.indx)
    loc.cox <- loc[indx.ppp[cox.indx],]

    write.csv(loc[indx.ppp[cox.indx],], file=paste(paste(output.loc, "loc.cox", sep="/"), count, "out", sep="."), row.names=FALSE)

    # Get the locations remaining after the homogeneous Poisson process
    X.q <- X[quad.indx,]
    n.q <- length(quad.indx)
    loc.q <- loc[quad.indx,]

    n.tot <- n.q + n.cox
    loc.tot <- rbind(loc.cox, loc.q)

    # Estimate the model by DWPR
    presence <- c(rep(1, n.cox), rep(0, n.q))
    #dd.tot <- deldir(loc.tot)
    #p.wt <- dd.tot$summary$dir.area
    p.wt <- rep(1/n.tot, n.tot)
    XX <- rbind(X.cox, X.q)
    dat <- data.frame(resp=round(presence/p.wt, 0), XX=XX, x=loc.tot$x, y=loc.tot$y, wt=p.wt)
    rownames(dat) <- 1:n.tot
    dat$id <- rownames(dat)

    dwpr <- glm(resp~XX, family=poisson(), weights=wt, data=dat)

    sink(file=paste(paste(output.loc, "dwpr", sep="/"), count, "out", sep="."), append=FALSE)
    print(dwpr$coefficients)
    sink()

    Z <- as.matrix(cbind(1, dat[,c('XX.1', 'XX.2', 'XX.3')]))


    #--------------------------------
    # Estimate the parameters of the Cox process regression model:
    spatial <- SRE(loc.tot)
    srecov <- SRE.covariance(sre=spatial, dwpr=dwpr, loc=loc.tot)
    Sp <- srecov$S %*% srecov$Lambda
    Sp <- sweep(Sp, 2, colMeans(Sp), '-')
    Sp2 <- sweep(Sp, 2, apply(Sp, 2, sd), '/')

    non.ortho.variational <- cox.variational(y=dat$resp, X=Z, wt=p.wt, S=Sp2, beta.start=dwpr$coefficients)
    non.ortho.variational.indep <- cox.variational.indep(y=dat$resp, X=Z, wt=p.wt, S=Sp2, beta.start=dwpr$coefficients)
    non.ortho.laplace <- cox.laplace(y=dat$resp, X=Z, wt=p.wt, S=Sp2, beta.start=dwpr$coefficients)

    sink(file=paste(paste(output.loc, "non.ortho", sep="/"), count, "out", sep="."), append=FALSE)
    cat(non.ortho)
    sink()


    #---------------------------------
    # Now estimate the model with random effects orthogonal to the fixed effects:
    srecov.o <- SRE.covariance(sre=spatial, dwpr=dwpr, loc=loc.tot, orthogonal=TRUE, X=as.data.frame(Z))
    Sp.o <- srecov.o$S %*% srecov.o$Lambda

    ortho.variational <- cox.variational(y=dat$resp, X=Z, wt=p.wt, S=Sp.o, beta.start=dwpr$coefficients)
    ortho.variational.indep <- cox.variational.indep(y=dat$resp, X=Z, wt=p.wt, S=Sp.o, beta.start=dwpr$coefficients)
    ortho.laplace <- cox.laplace(y=dat$resp, X=Z, wt=p.wt, S=Sp.o, beta.start=dwpr$coefficients)

    sink(file=paste(paste(output.loc, "ortho", sep="/"), count, "out", sep="."), append=FALSE)
    cat(ortho)
    sink()
  }
}
