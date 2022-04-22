// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sparse_plsa_3d
Rcpp::List sparse_plsa_3d(const Eigen::VectorXi& aye, const Eigen::VectorXi& jay, const Eigen::VectorXi& kay, const Eigen::VectorXd& qty, const int num_a, const int num_b, const int num_c, const int nz, const int niter);
RcppExport SEXP _HOPLSA_sparse_plsa_3d(SEXP ayeSEXP, SEXP jaySEXP, SEXP kaySEXP, SEXP qtySEXP, SEXP num_aSEXP, SEXP num_bSEXP, SEXP num_cSEXP, SEXP nzSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXi& >::type aye(ayeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi& >::type jay(jaySEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi& >::type kay(kaySEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type qty(qtySEXP);
    Rcpp::traits::input_parameter< const int >::type num_a(num_aSEXP);
    Rcpp::traits::input_parameter< const int >::type num_b(num_bSEXP);
    Rcpp::traits::input_parameter< const int >::type num_c(num_cSEXP);
    Rcpp::traits::input_parameter< const int >::type nz(nzSEXP);
    Rcpp::traits::input_parameter< const int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_plsa_3d(aye, jay, kay, qty, num_a, num_b, num_c, nz, niter));
    return rcpp_result_gen;
END_RCPP
}
// sparse_plsa_2d
Rcpp::List sparse_plsa_2d(const Eigen::VectorXi& aye, const Eigen::VectorXi& jay, const Eigen::VectorXd& qty, const int num_a, const int num_b, const int nz, const int niter);
RcppExport SEXP _HOPLSA_sparse_plsa_2d(SEXP ayeSEXP, SEXP jaySEXP, SEXP qtySEXP, SEXP num_aSEXP, SEXP num_bSEXP, SEXP nzSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXi& >::type aye(ayeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi& >::type jay(jaySEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type qty(qtySEXP);
    Rcpp::traits::input_parameter< const int >::type num_a(num_aSEXP);
    Rcpp::traits::input_parameter< const int >::type num_b(num_bSEXP);
    Rcpp::traits::input_parameter< const int >::type nz(nzSEXP);
    Rcpp::traits::input_parameter< const int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_plsa_2d(aye, jay, qty, num_a, num_b, nz, niter));
    return rcpp_result_gen;
END_RCPP
}
// sparse_plsa_4d
Rcpp::List sparse_plsa_4d(const Eigen::VectorXi& aye, const Eigen::VectorXi& jay, const Eigen::VectorXi& kay, const Eigen::VectorXi& ell, const Eigen::VectorXd& qty, const int num_a, const int num_b, const int num_c, const int num_d, const int nz, const int niter);
RcppExport SEXP _HOPLSA_sparse_plsa_4d(SEXP ayeSEXP, SEXP jaySEXP, SEXP kaySEXP, SEXP ellSEXP, SEXP qtySEXP, SEXP num_aSEXP, SEXP num_bSEXP, SEXP num_cSEXP, SEXP num_dSEXP, SEXP nzSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXi& >::type aye(ayeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi& >::type jay(jaySEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi& >::type kay(kaySEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi& >::type ell(ellSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type qty(qtySEXP);
    Rcpp::traits::input_parameter< const int >::type num_a(num_aSEXP);
    Rcpp::traits::input_parameter< const int >::type num_b(num_bSEXP);
    Rcpp::traits::input_parameter< const int >::type num_c(num_cSEXP);
    Rcpp::traits::input_parameter< const int >::type num_d(num_dSEXP);
    Rcpp::traits::input_parameter< const int >::type nz(nzSEXP);
    Rcpp::traits::input_parameter< const int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_plsa_4d(aye, jay, kay, ell, qty, num_a, num_b, num_c, num_d, nz, niter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HOPLSA_sparse_plsa_3d", (DL_FUNC) &_HOPLSA_sparse_plsa_3d, 9},
    {"_HOPLSA_sparse_plsa_2d", (DL_FUNC) &_HOPLSA_sparse_plsa_2d, 7},
    {"_HOPLSA_sparse_plsa_4d", (DL_FUNC) &_HOPLSA_sparse_plsa_4d, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_HOPLSA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
