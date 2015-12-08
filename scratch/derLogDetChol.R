# Rate of change of log determinant of Cholesky factor, w.r.t the original matrix (before Cholesky)

r <- nrow(P)

L <-t(chol(P))
Li <- solve(L)
dLi <- diag(Li)


D <- matrix(0, r, r)

for (j in 1:r) for (k in j:r) {
  der_diag_chol <- rep(0, r)

  for(l in j:r)
    der_diag_chol[l] <- diag(L)[l] * Li[l,j]*Li[l,k]/2

  D[j,k] <- sum(dLi * der_diag_chol)
}

D[lower.tri(D)] <- D[upper.tri(D)]



#----------------------------------
# Now want derivative of log determinant of Cholesky factor w.r.t. parameters beta and tau that appear in the original matrix:
# P <- tau*I_r + t(S) %*% diag(mu*wt*y) %*% S:


DP.Dbeta <- array(data=0, dim=c(r,r,p))
for (j in 1:p) {
  DP.Dbeta[,,j] <- t(S) %*% diag(mu*wt*X[,j]) %*% S
}

D2 <- rep(0, p)
for (j in 1:p) {
  D2[j] <- sum(D * DP.dbeta[,,j])
}


Dtau <- sum(diag(D))
