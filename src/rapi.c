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
    unsigned int *colsToPermutePtr = (unsigned int *)INTEGER(colsToPermute);
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

// try openmp backend with a simple parallel for loop in C with reduction(+:sum)
// #pragma omp parallel for reduction(+:sum)
SEXP R_Sum(SEXP x, SEXP nthreads)
{
    double *xptr = REAL(x);
    size_t n = Rf_length(x);
    int nthreads0 = INTEGER(nthreads)[0];
    double sum = 0;
    omp_set_num_threads(nthreads0);
    // int nthreads1 = omp_get_max_threads();
#pragma omp parallel for reduction(+:sum)
    for (size_t i = 0; i < n; i++)
    {
        sum += xptr[i];
        printf("thread id: %d/%d, i: %ld, sum: %f\n", omp_get_thread_num(),omp_get_num_threads(), i, sum);
    }
    // Rprintf("number of thread allocated: %d, requested: %d\n", nthreads1, nthreads0);
    return Rf_ScalarReal(sum);
}


// matrix add
SEXP R_MatrixAdd(SEXP x, SEXP y)
{
    size_t n = Rf_length(y);
    int nrow = Rf_nrows(x);
    int ncol = Rf_ncols(x);
    printf("nrow: %d, ncol: %d\n", nrow, ncol);
    SEXP rst = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
    double *xptr = REAL(x);
    double *yptr = REAL(y);
    double *rstptr = REAL(rst);
    omp_set_num_threads(4);

// #pragma omp parallel
// {
//     printf("thread %d/%d\n", omp_get_thread_num(), omp_get_max_threads());
// }

#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        rstptr[i] = xptr[i] + yptr[i];
        printf("thread id: %d, i: %ld, sum: %f\n", omp_get_thread_num(), i, rstptr[i]);
    }
    UNPROTECT(1);

    // #pragma omp parallel num_threads(2) // Set 2 threads for this parallel region
    // {
    //     int thread_id = omp_get_thread_num();
    //     // printf("Thread ID in parallel region with 2 threads: %d\n", thread_id);
    // }

    Rprintf("number of thread requested/allocated: %d/%d\n", 4L, omp_get_max_threads());
    return rst;
}

// return a list of vectors to R

SEXP R_ListOfVectors()
{
    SEXP list = PROTECT(Rf_allocVector(VECSXP, 3));

    double v1ptr[3];
    int v2ptr[4];
    char **v3ptr;
    v3ptr = (char **)malloc(5 * sizeof(char *));
    for (size_t i = 0; i < 5; i++)
    {
        v3ptr[i] = (char *)malloc(10 * sizeof(char));
    }   

    for(size_t i = 0; i < 3; i++)
    {
        v1ptr[i] = i * 1.0f;
    }

    for(size_t i = 0; i < 4; i++)
    {
        v2ptr[i] = i;
    }

    for(size_t i = 0; i < 5; i++)
    {
        sprintf(v3ptr[i], "string %ld", i);
    }

    SET_VECTOR_ELT(list, 0, DoublePtrToRVector(v1ptr, 3));
    SET_VECTOR_ELT(list, 1, IntPtrToRVector(v2ptr, 4));
    SET_VECTOR_ELT(list, 2, CharPtrToRVector(v3ptr, 5));

    UNPROTECT(1);

    return list;
}