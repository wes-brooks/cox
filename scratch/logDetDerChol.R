# Rate of change of log determinant of Cholesky factor, w.r.t the original matrix (before Cholesky)

r <- nrow(P)

L <-t(chol(P))
Li <- solve(L)
dLi <- diag(Li)


D <- matrix(0, r, r)

for (j in 1:r) for (k in j:r) {
  log_det_der_chol <- rep(0, r)

  for(l in j:r)
    log_det_der_chol[l] <- diag(L)[l] * Li[l,j]*Li[l,k]/2

  D[j,k] <- sum(dLi * log_det_der_chol)
}

D[lower.tri(D)] <- D[upper.tri(D)]
