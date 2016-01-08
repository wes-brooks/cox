#--------------------------
# Gradient of the lower likelihood bound w.r.t. V:
score.diagV <- function(logV, M, ltau, y, X, S, beta, wt) {
  #grad <- as.vector(colSums(sweep(S^2, 1, wt * mu * v / 2, '*')))
  eta <- as.vector(X %*% beta + S %*% M)
  mu <- exp(eta)
  tau <- exp(ltau)
  diagV <- exp(logV)

  v <- exp(VariationalVarIndep(diagV, S) / 2)

  grad.logV <- -as.vector(t(S^2) %*% as.matrix(wt * mu * v))*diagV / 2
  #grad.logV <- grad.logV + 1
  grad.logV <- grad.logV - tau*diagV / 2 # D_V
  #grad.logV <- grad.logV + DerLogDetCholIndep(diagV)
  grad.logV <- grad.logV + 1/2

  grad.logV
}



#--------------------------
# Gradient of the lower likelihood bound w.r.t. V:
score.variational <- function(M, logV, ltau, y, X, S, beta, wt) {
  log.diagV <- logV

  #grad <- as.vector(colSums(sweep(S^2, 1, wt * mu * v / 2, '*')))
  eta <- as.vector(X %*% beta + S %*% M)
  mu <- exp(eta)
  tau <- exp(ltau)
  diagV <- exp(log.diagV)

  v <- exp(VariationalVarIndep(diagV, S) / 2)

  # Gradient w.r.t. M:
  grad.M <- as.vector(t(S) %*% Diagonal(x=wt) %*% (y - mu*v)) - tau*M

  grad.logV <- -as.vector(t(S^2) %*% as.matrix(wt * mu * v))*diagV / 2
  #grad.logV <- grad.logV + 1
  grad.logV <- grad.logV - tau*diagV / 2 # D_V
  #grad.logV <- grad.logV + DerLogDetCholIndep(diagV)
  grad.logV <- grad.logV + 1/2

  c(grad.M, grad.logV)
}


# Gradient of the lower likelihood bound w.r.t. V:
score.fin.indep <- function(par, y, X, S, wt) {
  p <- ncol(X)
  r <- ncol(S)

  beta <- par[1:p]
  M <- par[(1:r)+p]
  log.diagV <- par[(1:r)+(r+p)]
  ltau <- tail(par, 1)

  #grad <- as.vector(colSums(sweep(S^2, 1, wt * mu * v / 2, '*')))
  eta <- as.vector(X %*% beta + S %*% M)
  mu <- exp(eta)
  tau <- exp(ltau)
  diagV <- exp(log.diagV)

  v <- exp(VariationalVarIndep(diagV, S) / 2)

  grad.logV <- -as.vector(t(S^2) %*% as.matrix(wt * mu * v))*diagV / 2
  #grad.logV <- grad.logV + 1
  grad.logV <- grad.logV - tau*diagV / 2
  grad.logV <- grad.logV + 1/2
  #grad.logV <- grad.logV + DerLogDetCholIndep(diagV)
  #grad.logV <- grad.logV / 2

  # Now gradient w.r.t. beta:
  grad.beta <- as.vector(t(X) %*% Diagonal(x=wt) %*% (y - mu*v))

  # Gradient w.r.t. M:
  grad.M <- as.vector(t(S) %*% Diagonal(x=wt) %*% (y - mu*v)) - tau*M

  # Gradient w.r.t. ltau:
  grad.ltau <- -tau/2 * (sum(M^2) + sum(diagV)) + r/2

  -c(grad.beta, grad.M, grad.logV, grad.ltau)
}
