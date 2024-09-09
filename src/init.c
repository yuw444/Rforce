#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* .call entry points */
extern SEXP R_Cummax(SEXP x);
extern SEXP R_ColsPermute(SEXP x, SEXP colsToPermute, SEXP seed);

static const R_CallMethodDef CallEntries[] = {
  {"R_Cummax", (DL_FUNC) &R_Cummax, 1},
  {"R_ColsPermute", (DL_FUNC) &R_ColsPermute, 3},
  {NULL, NULL, 0}
};

void R_init_RC(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
