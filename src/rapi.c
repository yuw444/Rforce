#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include "utils.h"
#include "convert.h"
#include <omp.h>

SEXP R_Cummax(SEXP x)
{
    double *xptr = REAL(x);
    size_t n = Rf_length(x);
    double *result = Cummax(xptr, n);
    SEXP res = DoublePtrToRVector(result, n);
    free(result);
    return res;
}

SEXP R_ColsPermute(SEXP x, SEXP colsToPermute, SEXP seed)
{

    int nrows = Rf_nrows(x);
    int ncols = Rf_ncols(x);
    // Rprintf("nrows: %d, ncols: %d\n", nrows, ncols);
    int ncolsToPermute = Rf_length(colsToPermute);
    double **xptr = RMatrixToDoublePtr(x);
    // for(int i = 0; i < nrows; i++)
    // {
    //     for(int j = 0; j < ncols; j++)
    //     {
    //         Rprintf("%f ", xptr[i][j]);
    //     }
    //     Rprintf("\n");
    // }
    int *colsToPermutePtr = INTEGER(colsToPermute);
    int seedInt = INTEGER(seed)[0];

    double **result = ColsPermute(xptr, nrows, ncols, colsToPermutePtr, ncolsToPermute, seedInt);
    // for(int i = 0; i < nrows; i++)
    // {
    //     for(int j = 0; j < ncols; j++)
    //     {
    //         Rprintf("%f ", result[i][j]);
    //     }
    //     Rprintf("\n");
    // }
    SEXP res = DoublePtrToRMatrix(result, nrows, ncols);
    return res;
}
