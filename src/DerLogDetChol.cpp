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
