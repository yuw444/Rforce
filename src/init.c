#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* .call entry points */
extern SEXP R_Rforce(
    SEXP splitFunctionIndex,
    SEXP interaction,
    SEXP designMatrixY,
    SEXP auxiliaryFeatures,
    SEXP varIDs,
    SEXP unitsOfCPIU,
    SEXP nTrees,
    SEXP maxDepth,
    SEXP minNodeSize,
    SEXP minGain,
    SEXP mtry,
    SEXP nsplits,
    SEXP seed);

extern SEXP R_ForestPredict(SEXP forestPtr, SEXP designMatrix);
extern SEXP R_SaveRforce(SEXP forestPtr, SEXP path);
extern SEXP R_LoadRforce(SEXP path);
extern SEXP R_PrintTree(SEXP forestPtr, SEXP treeIndex, SEXP filename);

static const R_CallMethodDef CallEntries[] = {
  {"R_Rforce", (DL_FUNC) &R_Rforce, 13},
  {"R_ForestPredict", (DL_FUNC) &R_ForestPredict, 2},
  {"R_SaveRforce", (DL_FUNC) &R_SaveRforce, 2},
  {"R_LoadRforce", (DL_FUNC) &R_LoadRforce, 1},
  {"R_PrintTree", (DL_FUNC) &R_PrintTree, 3},
  {NULL, NULL, 0}
};

void R_init_RC(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
