// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// perform_abc_rcpp_par
Rcpp::NumericMatrix perform_abc_rcpp_par(int num_particles, int num_iterations, double crown_age, double min_lin, double max_lin, std::vector<double> lambdas, double s, double obs_gamma, double obs_colless, double obs_num_lin);
RcppExport SEXP _pbdABC_perform_abc_rcpp_par(SEXP num_particlesSEXP, SEXP num_iterationsSEXP, SEXP crown_ageSEXP, SEXP min_linSEXP, SEXP max_linSEXP, SEXP lambdasSEXP, SEXP sSEXP, SEXP obs_gammaSEXP, SEXP obs_collessSEXP, SEXP obs_num_linSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type num_particles(num_particlesSEXP);
    Rcpp::traits::input_parameter< int >::type num_iterations(num_iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type crown_age(crown_ageSEXP);
    Rcpp::traits::input_parameter< double >::type min_lin(min_linSEXP);
    Rcpp::traits::input_parameter< double >::type max_lin(max_linSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type obs_gamma(obs_gammaSEXP);
    Rcpp::traits::input_parameter< double >::type obs_colless(obs_collessSEXP);
    Rcpp::traits::input_parameter< double >::type obs_num_lin(obs_num_linSEXP);
    rcpp_result_gen = Rcpp::wrap(perform_abc_rcpp_par(num_particles, num_iterations, crown_age, min_lin, max_lin, lambdas, s, obs_gamma, obs_colless, obs_num_lin));
    return rcpp_result_gen;
END_RCPP
}
// perform_abc_rcpp
Rcpp::NumericMatrix perform_abc_rcpp(int num_particles, int num_iterations, double crown_age, double min_lin, double max_lin, std::vector<double> lambdas, double s, double obs_gamma, double obs_colless, double obs_num_lin);
RcppExport SEXP _pbdABC_perform_abc_rcpp(SEXP num_particlesSEXP, SEXP num_iterationsSEXP, SEXP crown_ageSEXP, SEXP min_linSEXP, SEXP max_linSEXP, SEXP lambdasSEXP, SEXP sSEXP, SEXP obs_gammaSEXP, SEXP obs_collessSEXP, SEXP obs_num_linSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type num_particles(num_particlesSEXP);
    Rcpp::traits::input_parameter< int >::type num_iterations(num_iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type crown_age(crown_ageSEXP);
    Rcpp::traits::input_parameter< double >::type min_lin(min_linSEXP);
    Rcpp::traits::input_parameter< double >::type max_lin(max_linSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type obs_gamma(obs_gammaSEXP);
    Rcpp::traits::input_parameter< double >::type obs_colless(obs_collessSEXP);
    Rcpp::traits::input_parameter< double >::type obs_num_lin(obs_num_linSEXP);
    rcpp_result_gen = Rcpp::wrap(perform_abc_rcpp(num_particles, num_iterations, crown_age, min_lin, max_lin, lambdas, s, obs_gamma, obs_colless, obs_num_lin));
    return rcpp_result_gen;
END_RCPP
}
// test_simulations
Rcpp::List test_simulations(double birth, double death, double crown_age);
RcppExport SEXP _pbdABC_test_simulations(SEXP birthSEXP, SEXP deathSEXP, SEXP crown_ageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type birth(birthSEXP);
    Rcpp::traits::input_parameter< double >::type death(deathSEXP);
    Rcpp::traits::input_parameter< double >::type crown_age(crown_ageSEXP);
    rcpp_result_gen = Rcpp::wrap(test_simulations(birth, death, crown_age));
    return rcpp_result_gen;
END_RCPP
}
// sim_pbd_cpp
Rcpp::List sim_pbd_cpp(double la0, double mu0, double la1, double mu1, double trans_rate, double max_t, double max_num_species, int num_tries);
RcppExport SEXP _pbdABC_sim_pbd_cpp(SEXP la0SEXP, SEXP mu0SEXP, SEXP la1SEXP, SEXP mu1SEXP, SEXP trans_rateSEXP, SEXP max_tSEXP, SEXP max_num_speciesSEXP, SEXP num_triesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type la0(la0SEXP);
    Rcpp::traits::input_parameter< double >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< double >::type la1(la1SEXP);
    Rcpp::traits::input_parameter< double >::type mu1(mu1SEXP);
    Rcpp::traits::input_parameter< double >::type trans_rate(trans_rateSEXP);
    Rcpp::traits::input_parameter< double >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< double >::type max_num_species(max_num_speciesSEXP);
    Rcpp::traits::input_parameter< int >::type num_tries(num_triesSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_pbd_cpp(la0, mu0, la1, mu1, trans_rate, max_t, max_num_species, num_tries));
    return rcpp_result_gen;
END_RCPP
}
// sim_pbd_conditional_cpp
Rcpp::List sim_pbd_conditional_cpp(double la0, double mu0, double la1, double mu1, double trans_rate, double max_t, double min_num_species, double max_num_species, int num_tries);
RcppExport SEXP _pbdABC_sim_pbd_conditional_cpp(SEXP la0SEXP, SEXP mu0SEXP, SEXP la1SEXP, SEXP mu1SEXP, SEXP trans_rateSEXP, SEXP max_tSEXP, SEXP min_num_speciesSEXP, SEXP max_num_speciesSEXP, SEXP num_triesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type la0(la0SEXP);
    Rcpp::traits::input_parameter< double >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< double >::type la1(la1SEXP);
    Rcpp::traits::input_parameter< double >::type mu1(mu1SEXP);
    Rcpp::traits::input_parameter< double >::type trans_rate(trans_rateSEXP);
    Rcpp::traits::input_parameter< double >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< double >::type min_num_species(min_num_speciesSEXP);
    Rcpp::traits::input_parameter< double >::type max_num_species(max_num_speciesSEXP);
    Rcpp::traits::input_parameter< int >::type num_tries(num_triesSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_pbd_conditional_cpp(la0, mu0, la1, mu1, trans_rate, max_t, min_num_species, max_num_species, num_tries));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pbdABC_perform_abc_rcpp_par", (DL_FUNC) &_pbdABC_perform_abc_rcpp_par, 10},
    {"_pbdABC_perform_abc_rcpp", (DL_FUNC) &_pbdABC_perform_abc_rcpp, 10},
    {"_pbdABC_test_simulations", (DL_FUNC) &_pbdABC_test_simulations, 3},
    {"_pbdABC_sim_pbd_cpp", (DL_FUNC) &_pbdABC_sim_pbd_cpp, 8},
    {"_pbdABC_sim_pbd_conditional_cpp", (DL_FUNC) &_pbdABC_sim_pbd_conditional_cpp, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_pbdABC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
