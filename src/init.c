#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

#include <signal.h>
#include <execinfo.h>
#include <unistd.h>   // for STDERR_FILENO

/* ===== Crash handler for segfaults, aborts, etc. ===== */

static void crash_handler(int sig) {
  void *array[64];
  size_t size;

  fprintf(stderr,
          "\n=== Rforce native code crash: signal %d (e.g., segfault) ===\n",
          sig);

  // Get backtrace
  size = backtrace(array, 64);
  backtrace_symbols_fd(array, size, STDERR_FILENO);

  // Exit immediately; continuing is unsafe
  _exit(1);
}

static void install_crash_handler(void) {
  signal(SIGSEGV, crash_handler);  // segmentation fault
  signal(SIGABRT, crash_handler);  // abort()
  // you can add SIGFPE, SIGILL, etc. if you want
}

/* ===== .Call entry points ===== */

extern SEXP R_Rforce(
    SEXP splitFunctionIndex,
    SEXP interaction,
    SEXP designMatrixY,
    SEXP auxiliaryFeatures,
    SEXP varIDs,
    SEXP nUniqueVars,
    SEXP unitsOfCPIU,
    SEXP nTrees,
    SEXP maxDepth,
    SEXP minNodeSize,
    SEXP minGain,
    SEXP mtry,
    SEXP nsplits,
    SEXP nThreads,
    SEXP seed);

extern SEXP R_ForestPredict(SEXP forestPtr, SEXP designMatrix);
extern SEXP R_SaveRforce(SEXP forestPtr, SEXP path);
extern SEXP R_LoadRforce(SEXP path);
extern SEXP R_PrintTree(SEXP forestPtr, SEXP treeIndex, SEXP filename);

static const R_CallMethodDef CallEntries[] = {
  {"R_Rforce",        (DL_FUNC) &R_Rforce,        15},
  {"R_ForestPredict", (DL_FUNC) &R_ForestPredict,  2},
  {"R_SaveRforce",    (DL_FUNC) &R_SaveRforce,     2},
  {"R_LoadRforce",    (DL_FUNC) &R_LoadRforce,     1},
  {"R_PrintTree",     (DL_FUNC) &R_PrintTree,      3},
  {NULL, NULL, 0}
};

void R_init_Rforce(DllInfo *dll) {
  install_crash_handler();  // install segfault handler

  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
