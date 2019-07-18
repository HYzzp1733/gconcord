// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// ccista_conv
SEXP ccista_conv(const Eigen::MatrixXd& S, const Eigen::MappedSparseMatrix<double>& X0_, const Eigen::MatrixXd& LambdaMat, double lambda2, double epstol, int maxitr, int steptype);
RcppExport SEXP _gconcord_ccista_conv(SEXP SSEXP, SEXP X0_SEXP, SEXP LambdaMatSEXP, SEXP lambda2SEXP, SEXP epstolSEXP, SEXP maxitrSEXP, SEXP steptypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double>& >::type X0_(X0_SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type LambdaMat(LambdaMatSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< double >::type epstol(epstolSEXP);
    Rcpp::traits::input_parameter< int >::type maxitr(maxitrSEXP);
    Rcpp::traits::input_parameter< int >::type steptype(steptypeSEXP);
    rcpp_result_gen = Rcpp::wrap(ccista_conv(S, X0_, LambdaMat, lambda2, epstol, maxitr, steptype));
    return rcpp_result_gen;
END_RCPP
}
// ccfista_conv
SEXP ccfista_conv(const Eigen::MatrixXd& S, const Eigen::MappedSparseMatrix<double>& X0_, const Eigen::MatrixXd& LambdaMat, double lambda2, double epstol, int maxitr, int steptype);
RcppExport SEXP _gconcord_ccfista_conv(SEXP SSEXP, SEXP X0_SEXP, SEXP LambdaMatSEXP, SEXP lambda2SEXP, SEXP epstolSEXP, SEXP maxitrSEXP, SEXP steptypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double>& >::type X0_(X0_SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type LambdaMat(LambdaMatSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< double >::type epstol(epstolSEXP);
    Rcpp::traits::input_parameter< int >::type maxitr(maxitrSEXP);
    Rcpp::traits::input_parameter< int >::type steptype(steptypeSEXP);
    rcpp_result_gen = Rcpp::wrap(ccfista_conv(S, X0_, LambdaMat, lambda2, epstol, maxitr, steptype));
    return rcpp_result_gen;
END_RCPP
}
// ccorig_conv
SEXP ccorig_conv(const Eigen::MatrixXd& S, const Eigen::MatrixXd& X0_, const Eigen::MatrixXd& LambdaMat, double lambda2, double epstol, int maxitr);
RcppExport SEXP _gconcord_ccorig_conv(SEXP SSEXP, SEXP X0_SEXP, SEXP LambdaMatSEXP, SEXP lambda2SEXP, SEXP epstolSEXP, SEXP maxitrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X0_(X0_SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type LambdaMat(LambdaMatSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< double >::type epstol(epstolSEXP);
    Rcpp::traits::input_parameter< int >::type maxitr(maxitrSEXP);
    rcpp_result_gen = Rcpp::wrap(ccorig_conv(S, X0_, LambdaMat, lambda2, epstol, maxitr));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gconcord_ccista_conv", (DL_FUNC) &_gconcord_ccista_conv, 7},
    {"_gconcord_ccfista_conv", (DL_FUNC) &_gconcord_ccfista_conv, 7},
    {"_gconcord_ccorig_conv", (DL_FUNC) &_gconcord_ccorig_conv, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_gconcord(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
