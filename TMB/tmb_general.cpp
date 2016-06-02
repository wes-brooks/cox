// Likelihood lower bound for the variational approximation
// allow the random effects a general covariance matrix
#include <TMB.hpp>

template<class Type>
  Type objective_function<Type>::operator() ()
{
  // define the data objects
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_MATRIX(S);
  DATA_VECTOR(wt);
  //DATA_VECTOR(M);
  //DATA_VECTOR(beta);
  //DATA_SCALAR(ltau);

  // define the parameter objects
  PARAMETER_MATRIX(V);
  PARAMETER_VECTOR(M);
  PARAMETER_VECTOR(beta);
  PARAMETER(ltau);

  // define variables
  int r = S.cols();
  int n = X.rows();
  int p = X.cols();
  vector<Type> eta(n);
  vector<Type> v(n);
  Type cum;

  // compute Cholesky decomposition of V
  matrix<Type> L = V.llt().matrixL();

  // compute variance part of the Gaussian MGF
  for (int i=0; i<n; i++) {
    v(i) = 0.0;
    for (int j=0; j<r; j++) {
      cum = 0.0;
      for (int k=j; k<r; k++)
        cum += S(i, k) * L(k, j);

      v(i) += pow(cum, 2);
    }
  }

  // compute log of determinant of the Cholesky factor
  Type logDetCholV = 0;
  for (int j=0; j<r; j++)
    logDetCholV += log(L(j, j));

  // compute the linear predictor
  for (int i=0; i<n; i++) {
    eta(i) = 0.0;

    for (int j=0; j<p; j++)
      eta(i) += X(i, j) * beta(j);

    for (int j=0; j<r; j++)
      eta(i) += S(i, j) * M(j);
  }

  // sum of squared random effects
  Type M2 = 0;
  for (int j=0; j<r; j++)
    M2 += M(j) * M(j);

  // Calculate the negative log-likelihood
  Type nll = -(r*ltau - exp(ltau) * (M2 + V.trace())) / 2.0; // Expectation of prior on random effects
  nll -= r/2.0 * (1.0 + log(2.0*PI)) + logDetCholV; // Add the entropy term

  // Poisson GLM contribution to log-likelihood
  for (int i=0; i<n; i++)
    nll -= wt(i) * (y(i)*eta(i) - exp(eta(i) + v(i)/2.0));

  return nll;
}

