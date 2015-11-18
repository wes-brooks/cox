#---------------
# Function to simulate GRFs:
GRF <- function(loc, sigma=1, tau=0.1, nugget=0, fixed=null, ...) {
  covariate.model <- RMexp(var=sigma^2, scale=tau)
  if (nugget!=0)
    covariate.model <- covariate.model + RMnugget(var=nugget)

  dat <- RFsimulate(covariate.model, x=loc$x, y=loc$y, ...)
  dat
}
