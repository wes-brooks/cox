#-----------------------------
# Estimate the covariance of the spatial random effects, from residuals of a DWPR model
#' @export
SRE.covariance <- function(loc, dwpr, X=NULL, orthogonal=FALSE, n.1=16, n.res=3, bw.scale=1.5, res.scale=3) {
  complete <- FALSE

  while (!complete) {
    tryCatch({
      sre <- SRE(loc, n.1, n.res, bw.scale, res.scale)

      # Extract some info from the provided structure:
      S <- sre$SRE
      bins <- sre$bins
      M <- ncol(bins)

      if (orthogonal) {
        if (is.null(X)) stop("Error: you asked for orthogonal spatial effects but did not specify the fixed design matrix.")
        if (!is.data.frame(X) & !is.matrix(X)) stop("Error: the fixed design matrix is not provided as a data.frame or matrix.")
        if (nrow(X)!=nrow(S)) stop("Error: the fixed effects do not have the same number of rows as the spatial effects.")

        # Project spatial random effect matrix S to be orthogonal to the fixed effects:
        qrX.Q <- qr.Q(qr(X))
        P.1 <- t(qrX.Q) %*% S
        S <- S - qrX.Q %*% P.1
      }

      # Compute the method of moments variance Sigma.hat:
      D <- Matrix(0, nrow(loc), ncol(bins))
      for (i in 1:M) D[bins[,i]==1,i] <- dwpr$residuals[bins[,i]==1]
      D2 <- Matrix::colSums(D) / colSums(bins)
      Sigma.hat <- D2 %*% t(D2)
      diag(Sigma.hat) <- apply(D, 2, function(z) sum(z^2)) / colSums(bins)

      # Estimate the matrix K
      qrs <- qr(t(bins) %*% S / colSums(bins))
      qrs.R.inv <- solve(qr.R(qrs))
      qrs.Q <- qr.Q(qrs)
      K <- qrs.R.inv %*% t(qrs.Q) %*% Sigma.hat %*% qrs.Q %*% t(qrs.R.inv)
      Lambda <- chol(K)

      complete = TRUE
    }, error = function(e) {}, finally = n.res <- n.res - 1)

  }

  list(S=S, Lambda=Lambda, n.res=n.res+1)
}
