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



# SEEMS TO WORK!?!?!?!
sum(diag(Phi(Li %*% zed %*% t(Li))))
