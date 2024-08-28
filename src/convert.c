#include "convert.h"

// Convert R matrix to double**
double **RMatrixToDoublePtr(SEXP mat) {
  if (!Rf_isReal(mat) || !Rf_isMatrix(mat)) {
    Rf_error("Input must be a numeric matrix.");
  }

  printf("I'm here\n");
  int nrow = INTEGER(Rf_getAttrib(mat, R_DimSymbol))[0];
  int ncol = INTEGER(Rf_getAttrib(mat, R_DimSymbol))[1];

  double *data = REAL(mat);
  double **matrix = (double **) R_alloc(ncol, sizeof(double *));

  for (int j = 0; j < ncol; j++) {
    matrix[j] = &data[j * nrow];
  }

  return matrix;
}

// Convert R matrix to int**
int **RMatrixToIntPtr(SEXP mat) {
  if (!Rf_isInteger(mat) || !Rf_isMatrix(mat)) {
    Rf_error("Input must be an integer matrix.");
  }

  int nrow = INTEGER(Rf_getAttrib(mat, R_DimSymbol))[0];
  int ncol = INTEGER(Rf_getAttrib(mat, R_DimSymbol))[1];

  int *data = INTEGER(mat);
  int **matrix = (int **) R_alloc(ncol, sizeof(int *));

  for (int j = 0; j < ncol; j++) {
    matrix[j] = &data[j * nrow];
  }

  return matrix;
}

// Convert double** to R matrix
SEXP DoublePtrToRMatrix(double **matrix, int nrow, int ncol) {
  SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
  double *data = REAL(mat);

  for (int j = 0; j < ncol; j++) {
    for (int i = 0; i < nrow; i++) {
      data[i + nrow * j] = matrix[j][i];  // Access in column-major order
    }
  }

  UNPROTECT(1);
  return mat;
}

// Convert int** to R matrix
SEXP IntPtrToRMatrix(int **matrix, int nrow, int ncol) {
  SEXP mat = PROTECT(Rf_allocMatrix(INTSXP, nrow, ncol));
  int *data = INTEGER(mat);

  for (int j = 0; j < ncol; j++) {
    for (int i = 0; i < nrow; i++) {
      data[i + nrow * j] = matrix[j][i];  // Access in column-major order
    }
  }

  UNPROTECT(1);
  return mat;
}

// Convert R vector to double*
double *RVectorToDoublePtr(SEXP vec) {
  if (!Rf_isReal(vec)) {
    Rf_error("Input must be a numeric vector.");
  }

  return REAL(vec);
}

// Convert R vector to int*
int *RVectorToIntPtr(SEXP vec) {
  if (!Rf_isInteger(vec)) {
    Rf_error("Input must be an integer vector.");
  }

  return INTEGER(vec);
}

// Convert double* to R vector
SEXP DoublePtrToRVector(double *vector, int length) {
  SEXP vec = PROTECT(Rf_allocVector(REALSXP, length));
  double *data = REAL(vec);

  for (int i = 0; i < length; i++) {
    data[i] = vector[i];
  }

  UNPROTECT(1);
  return vec;
}

// Convert int* to R vector
SEXP IntPtrToRVector(int *vector, int length) {
  SEXP vec = PROTECT(Rf_allocVector(INTSXP, length));
  int *data = INTEGER(vec);

  for (int i = 0; i < length; i++) {
    data[i] = vector[i];
  }

  UNPROTECT(1);
  return vec;
}
