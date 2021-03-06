// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// DerLogDetChol
Eigen::MatrixXd DerLogDetChol(const Eigen::MatrixXd L);
RcppExport SEXP cox_DerLogDetChol(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type L(LSEXP);
    __result = Rcpp::wrap(DerLogDetChol(L));
    return __result;
END_RCPP
}
// DerLogDetCholIndep
Eigen::VectorXd DerLogDetCholIndep(const Eigen::VectorXd diagV);
RcppExport SEXP cox_DerLogDetCholIndep(SEXP diagVSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type diagV(diagVSEXP);
    __result = Rcpp::wrap(DerLogDetCholIndep(diagV));
    return __result;
END_RCPP
}
// DerLogDetCholTau
double DerLogDetCholTau(const Eigen::MatrixXd L);
RcppExport SEXP cox_DerLogDetCholTau(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type L(LSEXP);
    __result = Rcpp::wrap(DerLogDetCholTau(L));
    return __result;
END_RCPP
}
// DerLogDetCholBeta
Eigen::VectorXd DerLogDetCholBeta(const Eigen::MatrixXd L, const Eigen::MatrixXd S, const Eigen::MatrixXd X, const Eigen::VectorXd beta, const Eigen::VectorXd u, const Eigen::VectorXd wt);
RcppExport SEXP cox_DerLogDetCholBeta(SEXP LSEXP, SEXP SSEXP, SEXP XSEXP, SEXP betaSEXP, SEXP uSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type L(LSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type S(SSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type wt(wtSEXP);
    __result = Rcpp::wrap(DerLogDetCholBeta(L, S, X, beta, u, wt));
    return __result;
END_RCPP
}
// VariationalVar
NumericVector VariationalVar(const Eigen::MatrixXd cholV, const Eigen::MatrixXd S);
RcppExport SEXP cox_VariationalVar(SEXP cholVSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type cholV(cholVSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type S(SSEXP);
    __result = Rcpp::wrap(VariationalVar(cholV, S));
    return __result;
END_RCPP
}
// VariationalVarIndep
NumericVector VariationalVarIndep(const Eigen::VectorXd diagV, const Eigen::MatrixXd S);
RcppExport SEXP cox_VariationalVarIndep(SEXP diagVSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type diagV(diagVSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type S(SSEXP);
    __result = Rcpp::wrap(VariationalVarIndep(diagV, S));
    return __result;
END_RCPP
}
// VariationalScore
Eigen::MatrixXd VariationalScore(const Eigen::VectorXd mu, const Eigen::VectorXd wt, double tau, const Eigen::VectorXd v, const Eigen::MatrixXd V, const Eigen::MatrixXd S);
RcppExport SEXP cox_VariationalScore(SEXP muSEXP, SEXP wtSEXP, SEXP tauSEXP, SEXP vSEXP, SEXP VSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type v(vSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type V(VSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type S(SSEXP);
    __result = Rcpp::wrap(VariationalScore(mu, wt, tau, v, V, S));
    return __result;
END_RCPP
}
// VariationalScoreLogV
Eigen::MatrixXd VariationalScoreLogV(const Eigen::VectorXd mu, const Eigen::VectorXd wt, double tau, const Eigen::VectorXd v, const Eigen::MatrixXd V, const Eigen::MatrixXd S);
RcppExport SEXP cox_VariationalScoreLogV(SEXP muSEXP, SEXP wtSEXP, SEXP tauSEXP, SEXP vSEXP, SEXP VSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type v(vSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type V(VSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type S(SSEXP);
    __result = Rcpp::wrap(VariationalScoreLogV(mu, wt, tau, v, V, S));
    return __result;
END_RCPP
}
