#include "convert.h"

/**
 * 1 2 3 4
 * 5 6 7 8
 * 9 10 11 12
 */

// Convert column-major matrix to row-major matrix
void ColMajorToRowMajor(double *colMajor, double *rowMajor, int nrow, int ncol)
{
  for (int j = 0; j < ncol; j++)
  {
    for (int i = 0; i < nrow; i++)
    {
      rowMajor[j + ncol * i] = colMajor[i + nrow * j]; // Access in row-major order
    }
  }
}

// Convert row-major matrix to column-major matrix
void RowMajorToColMajor(double *rowMajor, double *colMajor, int nrow, int ncol)
{
  for (int j = 0; j < ncol; j++)
  {
    for (int i = 0; i < nrow; i++)
    {
      colMajor[i + nrow * j] = rowMajor[j + ncol * i]; // Access in column-major order
    }
  }
}

// Convert R matrix to double**
double **RMatrixToDoublePtr(SEXP mat)
{

  if (!Rf_isReal(mat) || !Rf_isMatrix(mat))
  {
    Rf_error("Input must be a numeric matrix.");
  }
  int nrow = Rf_nrows(mat);
  int ncol = Rf_ncols(mat);

  double *data = REAL(mat);
  double **matrix = Allocate2DArray(nrow, ncol);
  for (int i = 0; i < nrow; i++)
  {
    for (int j = 0; j < ncol; j++)
    {
      matrix[i][j] = ISNA(data[i + j * nrow]) ? NA_DOUBLE : data[i + j * nrow];
    }
    // if(i < 5)
    //   PrintArrayDouble(matrix[i], ncol);
  }

  return matrix;
}

// Convert R matrix to int**
int **RMatrixToIntPtr(SEXP mat)
{
  if (!Rf_isInteger(mat) || !Rf_isMatrix(mat))
  {
    Rf_error("Input must be an integer matrix.");
  }

  int nrow = Rf_nrows(mat);
  int ncol = Rf_ncols(mat);

  int *data = INTEGER(mat);
  int **matrix = Allocate2DArrayInt(nrow, ncol);
  for (int i = 0; i < nrow; i++)
  {
    for (int j = 0; j < ncol; j++)
    {
      matrix[i][j] = data[i + j * nrow];
    }
  }

  return matrix;
}

// Convert double** to R matrix
SEXP DoublePtrToRMatrix(double **matrix, int nrow, int ncol)
{
  SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
  double *data = REAL(mat);

  for (int j = 0; j < ncol; j++)
  {
    for (int i = 0; i < nrow; i++)
    {
      data[i + nrow * j] = matrix[i][j]; // Access in column-major order
    }
  }

  UNPROTECT(1);
  return mat;
}

// Convert int** to R matrix
SEXP IntPtrToRMatrix(int **matrix, int nrow, int ncol)
{
  SEXP mat = PROTECT(Rf_allocMatrix(INTSXP, nrow, ncol));
  int *data = INTEGER(mat);

  for (int j = 0; j < ncol; j++)
  {
    for (int i = 0; i < nrow; i++)
    {
      data[i + nrow * j] = matrix[i][j]; // Access in column-major order
    }
  }

  UNPROTECT(1);
  return mat;
}

// Convert R vector to double*
double *RVectorToDoublePtr(SEXP vec)
{
  if (!Rf_isReal(vec))
  {
    Rf_error("Input must be a numeric vector.");
  }

  return REAL(vec);
}

// Convert R vector to int*
int *RVectorToIntPtr(SEXP vec)
{
  if (!Rf_isInteger(vec))
  {
    Rf_error("Input must be an integer vector.");
  }

  return INTEGER(vec);
}

// Convert double* to R vector
SEXP DoublePtrToRVector(double *vector, int length)
{
  SEXP vec = PROTECT(Rf_allocVector(REALSXP, length));
  double *data = REAL(vec);

  for (int i = 0; i < length; i++)
  {
    data[i] = vector[i];
  }

  UNPROTECT(1);
  return vec;
}

// Convert int* to R vector
SEXP IntPtrToRVector(int *vector, int length)
{
  SEXP vec = PROTECT(Rf_allocVector(INTSXP, length));
  int *data = INTEGER(vec);

  for (int i = 0; i < length; i++)
  {
    data[i] = vector[i];
  }

  UNPROTECT(1);
  return vec;
}

// convert char** to R character vector
SEXP CharPtrToRVector(char **vector, int length)
{
  SEXP vec = PROTECT(Rf_allocVector(STRSXP, length));

  for (int i = 0; i < length; i++)
  {
    SET_STRING_ELT(vec, i, Rf_mkChar(vector[i]));
  }

  UNPROTECT(1);
  return vec;
}
