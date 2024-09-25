#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include "utils.h"

// Convert R matrix to double**
double **RMatrixToDoublePtr(SEXP mat);

// Convert R matrix to int**
int **RMatrixToIntPtr(SEXP mat);

// Convert R vector to double*
double *RVectorToDoublePtr(SEXP vec);

// Convert R vector to int*
int *RVectorToIntPtr(SEXP vec);

// Convert double** to R matrix
SEXP DoublePtrToRMatrix(double **matrix, int nrow, int ncol);

// Convert int** to R matrix
SEXP IntPtrToRMatrix(int **matrix, int nrow, int ncol);

// Convert double* to R vector
SEXP DoublePtrToRVector(double *vector, int length);

// Convert int* to R vector
SEXP IntPtrToRVector(int *vector, int length);

// Convert char** to R vector
SEXP CharPtrToRVector(char **vector, int length);
