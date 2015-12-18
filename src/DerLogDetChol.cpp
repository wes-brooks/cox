#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]

// This function computes the derivative of the log determinant of a Cholesky factor
// with respect to regression coefficients beta and precision component tau

// [[Rcpp::export]]
Eigen::MatrixXd DerLogDetChol(const Eigen::MatrixXd L)
{
  int r = L.cols();

  Eigen::MatrixXd Li = L.inverse();
  Eigen::VectorXd dLi = Li.diagonal();
  Eigen::MatrixXd D(r, r);
  D.setConstant(0);

  for (int j=0; j<r; j++)
  {
    for (int k=j; k<r; k++)
    {
      for(int l=j; l<r; l++)
      {
        if (j!=k)
          D(j,k) += dLi(l) * L.diagonal()(l) * Li(l,j) * Li(l,k);
        else
          D(j,k) += dLi(l) * L.diagonal()(l) * Li(l,j) * Li(l,k) / 2;
      }
    }
  }

  return D;
}



// [[Rcpp::export]]
Eigen::VectorXd DerLogDetCholIndep(const Eigen::VectorXd diagV)
{
  int r = diagV.count();

  Eigen::VectorXd D(r);
  D.setConstant(0);

  for (int j=0; j<r; j++)
  {
    D(j) = 1 / diagV(j) / 2;
  }

  return D;
}




// [[Rcpp::export]]
Eigen::VectorXd DerLogDetChol3(const Eigen::MatrixXd L)
{
  int r = L.cols();
  int size = r * (r+1) / 2;
  Eigen::MatrixXd Li = L.inverse();
  Eigen::VectorXd D(size);
  D.setConstant(0);

  int indx = 0;
  for (int j=0; j<r; j++)
  {
    for (int k=j; k<r; k++)
    {
      D(indx) += pow(Li(k, j), 2) / 2;
    }
    indx += 1;

    for (int i=j+1; i<r; i++)
    {
      for (int k=i; k<r; k++)
        D(indx) += Li(k, i) * Li(k, j);

      indx += 1;
    }
  }

  return(D);
}



// [[Rcpp::export]]
double DerLogDetCholTau(const Eigen::MatrixXd L)
{
  int r = L.cols();
  Eigen::MatrixXd Li = L.inverse();
  double result=0;

  int indx = 0;
  for (int j=0; j<r; j++)
  {
    for (int k=j; k<r; k++)
      result += pow(Li(k, j), 2) / 2;
  }

  return(result);
}

//
//
// DerLogDetCholBeta <- function(P, S, X, beta, u, wt) {
//   L <- t(chol(P))
//   Li <- solve(L)
//   r <- ncol(P)
//   p <- ncol(X)
//   n <- nrow(X)
//   D <- matrix(0, p, r)
//
//   mu <- exp(X %*% beta + S %*% u)
//
//
//   for (j in 1:r) {
//     AS <- vector()
//     for (i in 1:n)
//       AS <- c(AS, wt[i] * mu[i] * sum(Li[j, 1:j] * S[i, 1:j])^2)
//
//
//       D[,j] <- L[j,j] * apply(X, 2, function(z) sum(z * AS) / 2)
//   }
//
//   apply(D, 1, function(z) sum(z / diag(L)))
// }


// [[Rcpp::export]]
Eigen::VectorXd DerLogDetCholBeta(const Eigen::MatrixXd L, const Eigen::MatrixXd S, const Eigen::MatrixXd X, const Eigen::VectorXd beta, const Eigen::VectorXd u, const Eigen::VectorXd wt)
{
  int r = L.cols();
  int p = X.cols();
  int n = X.rows();

  Eigen::MatrixXd Li = L.inverse();
  Eigen::MatrixXd D(p, r);
  D.setConstant(0);

  Eigen::VectorXd nu(n);
  nu = X*beta + S*u;

  Eigen::VectorXd result(p);

  for (int j=0; j<r; j++)
  {
    Eigen::VectorXd LiS(n);
    LiS.setConstant(0);

    for (int i=0; i<n; i++)
    {
      for (int k=0; k<=j; k++)
        LiS(i) += Li(j, k) * S(i, k);

      LiS(i) = pow(LiS(i), 2) * wt(i) * exp(nu(i));
    }

    for (int k=0; k<p; k++)
      D(k, j) = (LiS.array() * X.col(k).array()).sum() / 2;
  }

  for (int j=0; j<p; j++)
    result(j) = D.row(j).sum();

  return(result);
}

//
// #----------------------------------
// # Now want derivative of log determinant of Cholesky factor w.r.t. parameters beta and tau that appear in the original matrix:
// # P <- tau*I_r + t(S) %*% diag(mu*wt) %*% S:
//
//
// DP.Dbeta <- array(data=0, dim=c(r,r,p))
// for (j in 1:p) {
//   DP.Dbeta[,,j] <- t(S) %*% diag(mu*wt*X[,j]) %*% S
// }
//
// D2 <- rep(0, p)
// for (j in 1:p) {
//   D2[j] <- sum(D * DP.dbeta[,,j])
// }
//
//
// Dtau <- sum(diag(D))
