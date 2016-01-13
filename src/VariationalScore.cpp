#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
NumericVector VariationalVar(const Eigen::MatrixXd cholV, const Eigen::MatrixXd S)
{
  int n = S.rows();
  NumericVector v(n);

  for (int i=0; i<n; i++)
    v(i) = (S.row(i) * cholV).array().pow(2).sum();

  return v;
}


// [[Rcpp::export]]
NumericVector VariationalVarIndep(const Eigen::VectorXd diagV, const Eigen::MatrixXd S)
{
  int n = S.rows();
  int p = S.cols();
  NumericVector v(n);

  for (int i=0; i<n; i++)
  {
    v(i) = 0;

    for (int j=0; j<p; j++)
      v(i) += pow(S(i, j), 2) * diagV(j);
  }

  return v;
}




// [[Rcpp::export]]
Eigen::MatrixXd VariationalScore(const Eigen::VectorXd mu, const Eigen::VectorXd wt, double tau, const Eigen::VectorXd v, const Eigen::MatrixXd V, const Eigen::MatrixXd S)
{
  int p = S.cols();

  // Gradient of the lower likelihood bound:
  Eigen::VectorXd diag = mu.array() * wt.array() * v.array();
  Eigen::MatrixXd grad = -S.transpose() * diag.asDiagonal() * S;

  grad += V.inverse();
  for (int i=0; i<p; i++)
    grad.diagonal()[i] -= tau;
  grad *= 0.5;

  return grad;
}



// [[Rcpp::export]]
Eigen::MatrixXd VariationalScoreLogV(const Eigen::VectorXd mu, const Eigen::VectorXd wt, double tau, const Eigen::VectorXd v, const Eigen::MatrixXd V, const Eigen::MatrixXd S)
{
  int p = S.cols();

  // Gradient of the lower likelihood bound:
  Eigen::VectorXd diag = mu.array() * wt.array() * v.array();
  Eigen::MatrixXd grad = -S.transpose() * diag.asDiagonal() * S;

  for (int i=0; i<p; i++)
    grad.diagonal()[i] -= tau;

  // This matrix has ones on the diagonal and twos off-diagonal.
  // The twos double the gradient off-diagonal because the variance-covariance matrix is symmetric
  // and we will only keep the upper triangle.
  Eigen::MatrixXd ones;
  ones.setOnes(p, p);
  ones = ones * 2;
  for (int i=0; i<p; i++)
    ones.diagonal()(i) = 1;

  grad = grad.array() * V.array() * ones.array();
  grad *= 0.5;

  return grad;
}

