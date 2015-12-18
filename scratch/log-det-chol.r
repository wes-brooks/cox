DerLogDetChol <- function(L) {
  r <- ncol(L)
  Li <- solve(L)
  res <- vector()

  for (j in 1:r) for (i in j:r) {
    x <- sum(tail(Li[,i], r-i+1) * tail(Li[,j], r-i+1) * tail(diag(L), r-i+1)) / 2

    res <- c(res, x)
  }

  res
}

Phi <- function(M) {
  res <- matrix(0, nrow(M), ncol(M))
  low.indx <- which(lower.tri(M))

  res[low.indx] <- M[low.indx]
  diag(res) <- diag(M) / 2

  res
}


DerLogDetChol2 <- function(P, zed) {
  A <- t(chol(P))
  r <- ncol(P)
  D <- matrix(0, r, r)

  for (i in 1:r) {
    indx <- which(1:r < i)
    jindx <- which(1:r > i)
    A[i,i] <- sqrt(P[i,i] - sum(A[i,indx]^2))
    D[i,i] <- (zed[i,i] - sum(2*D[i,indx]*A[i,indx])) / A[i,i] / 2

    for (j in jindx) {
      A[j,i] <- (P[j,i] - sum(A[j,indx] * A[i,indx])) / A[i,i]
      tmp <- zed[j,i] - sum(D[j,indx] * A[i,indx] + A[j,indx] * D[i, indx])
      D[j,i] <- tmp / A[i,i] - A[j,i] * (D[i,i] / A[i,i])
    }
  }

  D
}


DerLogDetChol4 <- function(P) {
  A <- t(chol(P))
  Ai <- solve(A)
  r <- ncol(P)
  D <- matrix(0, r, r)

  for (j in 1:r) {
    D[j,j] <- sum(Ai[j:r, j]^2) / 2

    indx <- which(1:r > j)
    for (i in indx) {
      D[i, j] <- sum(Ai[i:r, i] * Ai[i:r, j])
    }
  }

  D
}

DerLogDetCholTau <- function(P) {
  L <- t(chol(P))
  Li <- solve(L)
  r <- ncol(P)
  D <- vector()

  for (j in 1:r) {
    D <- c(D, L[j,j] * sum(Li[j:r,j]^2) / 2)
  }

  sum(D / diag(L))
}


DerLogDetCholTau2 <- function(P) {
  L <- t(chol(P))
  Li <- solve(L)
  r <- ncol(P)
  D <- vector()

  for (j in 1:r) {
    D <- c(D, sum(Li[j:r,j]^2) / 2)
  }

  sum(D)
}



DerLogDetCholBeta2 <- function(l, S, X, beta, u, wt) {
  Li <- solve(L)
  r <- ncol(P)
  p <- ncol(X)
  n <- nrow(X)
  D <- matrix(0, p, r)

  mu <- exp(X %*% beta + S %*% u)


  for (j in 1:r) {
    AS <- vector()
    for (i in 1:n)
      AS <- c(AS, wt[i] * mu[i] * sum(Li[j, 1:j] * S[i, 1:j])^2)


    D[,j] <- apply(X, 2, function(z) sum(z * AS) / 2)
  }

  rowSums(D)
}




data(mtcars)
X <- cbind(1, as.matrix(mtcars[2:11]))
y <- mtcars$mpg
n <- nrow(X)
S <- cbind(rnorm(n), rnorm(n), rnorm(n), rnorm(n))
r <- ncol(S)
u <- rnorm(r)
wt <- runif(n)
tau <- 6

beta <- lm(mpg~., data=mtcars, weights=wt)$coefficients

mu <- as.vector(exp(X%*%beta + S%*%u))

P <- diag(x=rep(tau, r)) + t(S) %*% diag(x=wt*mu) %*% S
L <- t(chol(P))
Li <- solve(L)

# Change beta and see change in log determinant of Cholesky P:
b2 <- beta
b2[1] <- b2[1] + 1/1000
m2 <- as.vector(exp(X %*% b2 + S %*% u))
P2 <- diag(x=rep(tau, r)) + t(S) %*% diag(x=wt*m2) %*% S
sum(log(diag(chol(P2)))) - sum(log(diag(chol(P))))


zed <- matrix(0, ncol(X), ncol(X))
zed[5,4] <- zed[4,5] <- 1

# SEEMS TO WORK!?!?!?!
sum(diag(Phi(Li %*% zed %*% t(Li))))

# COMPARE TO:
sum(log(diag(chol(P + zed)))) - sum(log(diag(chol(P))))
