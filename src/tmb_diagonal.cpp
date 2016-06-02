// Likelihood lower bound for the variational approximation
// the random effects have a diagonal covariance matrix
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
  PARAMETER_VECTOR(logV);
  PARAMETER_VECTOR(M);
  PARAMETER_VECTOR(beta);
  PARAMETER(ltau);

  // define variables
  int r = S.cols();
  int n = X.rows();
  int p = X.cols();
  vector<Type> eta(n);
  vector<Type> v(n);

  // calculate the variance part of the Gaussian MGF
  for (int i=0; i<n; i++) {
    v(i) = 0.0;

    for (int j=0; j<r; j++)
      v(i) += S(i, j) * S(i, j) * exp(logV(j));
  }

  // calculate the linear predictor
  for (int i=0; i<n; i++) {
    eta(i) = 0.0;

    for (int j=0; j<p; j++)
      eta(i) += X(i, j) * beta(j);

    for (int k=0; k<r; k++)
      eta(i) += S(i, k) * M(k);
  }

  // sum of squared random effects
  Type M2 = 0;
  for (int j=0; j<r; j++)
    M2 += M(j) * M(j);

  // negative log-lik of the not-Poisson part:
  Type nll = -(r*ltau - exp(ltau)*(M2 + exp(logV).sum())) / 2.0; // Add expectation of the log prior on the random effects
  nll -= r/2.0 * (1.0 + log(2.0*PI)) + logV.sum() / 2.0; // Add the entropy term

  // Poisson part of the likelihood
  for (int i=0; i<n; i++)
    nll -= wt(i) * (y(i)*eta(i) - exp(eta(i) + v(i)/2.0));

  return nll;
}
