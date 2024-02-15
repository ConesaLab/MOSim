// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <random>
#include <iostream>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector random_unif_interval(int size, int max_val){
  
  // vector to store the generated random number
  Rcpp::NumericVector index(size);
  
  // initialize the random generator
  std::random_device seed;
  std::mt19937 rng(seed());
  
  // uncomment this line if std::random_device does not work properly
  //std::default_random_engine rng(time(0));
  
  std::uniform_int_distribution<> unif_c(1, max_val);
  for(int i=0; i<size; ++i){
    index[i] = unif_c(rng);
  }
  return(index);
  
}

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// random_unif_interval
NumericVector random_unif_interval(int size, int max_val);
RcppExport SEXP _sparsim_random_unif_interval(SEXP sizeSEXP, SEXP max_valSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
  Rcpp::traits::input_parameter< int >::type max_val(max_valSEXP);
  rcpp_result_gen = Rcpp::wrap(random_unif_interval(size, max_val));
  return rcpp_result_gen;
  END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
  {"_sparsim_random_unif_interval", (DL_FUNC) &_sparsim_random_unif_interval, 2},
  {NULL, NULL}
};

RcppExport void R_init_SPARSim(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
