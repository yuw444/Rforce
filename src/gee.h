#ifndef GEE_H
#define GEE_H

#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <stdlib.h>
#include <setjmp.h>
#include <dirent.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <stdarg.h>
#include "utils.h"

#define MAX_ROWS 100
#define PERMANENT 1
#define EPHEMERAL 0
#define MAX_COVLAG 30
#define NO_ERROR 0
#define UNKNOWN_FAILURE 1
#define NO_MEM_MATSTRUCT 2
#define NO_MEM_MATDATA 3
#define SPLIT_FAIL 4
#define MATMULT_NONCONFORMITY 5
#define MATADD_NONCONFORMITY 6
#define CCHOL_FAIL 7
#define CORNER_FAIL 8
#define EXCEED_MAX_COVLAG 9
#define BAD_TOEPLITZ_ARG 10
#define PLUG_FAIL 11
#define PX1XPXQ_ARG1_BAD 12
#define PX1XPXQ_CONFORMITY 13
#define PXQDPX1_ARG1_BAD 14
#define PXQDPX1_CONFORMITY 15
#define CCHOL_NOT_SQUARE 16
#define MATREAD_OPEN_FAIL 17
#define MATREAD_NOT_RECTANGLE 18
#define INV_FAILURE 19

extern int _GSL_ERROR_FLAG;

// typedef struct matrix
// {
//     int nrows, ncols;
//     double *data;
//     int permanence;
// } MATRIX;
// /* element reference is handled principally by MEL */
// #define ELREF(matp, s1, s2) ((matp)->data) + (s2) + ((s1) * (matp->ncols))
// #define MEL(X, i, j) (*(ELREF((X), (i), (j))))
// #define get_nelem(x) (((x)->nrows) * ((x)->ncols))

typedef struct matrix
{
    int nrows, ncols;
    double **data;
    int permanence;
} MATRIX;
/* element reference is handled principally by MEL */
#define get_nelem(x) (((x)->nrows) * ((x)->ncols))

#define is_permanent(x) (x)->permanence == PERMANENT
#define is_ephemeral(x) (x)->permanence == EPHEMERAL
#define make_permanent(x) (x)->permanence = PERMANENT;
#define make_ephemeral(x) (x)->permanence = EPHEMERAL;

void VC_GEE_destroy_matrix(MATRIX **mat);

void VC_GEE_destroy_double_matrix(MATRIX ***mat, int nmat);

void VC_GEE_destroy_matrices(int nmat, ...);

void free_if_ephemeral(MATRIX **mat);

/*The Sdblptr is concated by column to form a long array*/
#define from_S(Sdblptr, Srowintptr, Scolintptr, Matptr)                           \
    Matptr = VC_GEE_create_matrix((int)*Srowintptr, (int)*Scolintptr, EPHEMERAL); \
    {                                                                             \
        int i, j, Scol, Srow;                                                     \
        double *Sload;                                                            \
        Scol = *Scolintptr;                                                       \
        Srow = *Srowintptr;                                                       \
        Sload = Sdblptr;                                                          \
        for (i = 0; i < Srow; i++)                                                \
        {                                                                         \
            for (j = 0; j < Scol; j++)                                            \
            {                                                                     \
                (Matptr->data)[i][j] = (double)*(Sload++);                           \
            }                                                                     \
        }                                                                         \
    }
/* end define |from_S| */

#define to_S(Matptr, Sdblptr)                   \
    {                                           \
        int i, j;                               \
        double *Sload;                          \
        Sload = Sdblptr;                        \
        for (i = 0; i < Matptr->nrows; i++)     \
        {                                       \
            for (j = 0; j < Matptr->ncols; j++) \
            {                                   \
                *(Sload++) = (Matptr->data)[i][j]; \
            }                                   \
        }                                       \
    }

MATRIX *VC_GEE_create_matrix(int, int, int),
    /*create a brand new copy of the matrix*/
    *VC_GEE_matcopy(MATRIX **),
    /*get multiple rows of the matrix*/
    *VC_GEE_extract_rows(MATRIX **, int, int),
    /*element-wise add two matrices, destroy the parameter matrices as needed, return a temporal new matrix*/
    *VC_GEE_matadd(MATRIX **, MATRIX **),
    /*element-wise substract two matrices, destroy the parameter matrices as needed, return a temporal new matrix*/
    *VC_GEE_matsub(MATRIX **, MATRIX **),
    /*dot product of two confront matrices, destroy the parameter matrices as needed, return a temporal new matrix*/
    *VC_GEE_matmult(MATRIX **, MATRIX **),
    /*the transpose of the matrix, destroy the parameter matrix as needed, return a temporal new matrix*/
    *VC_GEE_transp(MATRIX **),
    /*return a temporal 1 column of 1s*/
    *VC_GEE_col_1s(int),
    /*elementwise absolute value of a matrix, destroy the parameter matrix as needed, return a temporal new matrix */
    *VC_GEE_matabs(MATRIX **),
    /*elementwise exp of a matrix, destroy the parameter matrix as needed, return a temporal new matrix */
    *VC_GEE_matexp(MATRIX **),
    /*multiple each row of matrix by the correspondent elment in the column vector, destroy the parameter matrix as needed, return a temporal new matrix*/
    *VC_GEE_px1_times_pxq(MATRIX **, MATRIX **),
    /*divide each column of matrix by the correspondent elment in the row vector, destroy the parameter matrix as needed, return a temporal new matrix*/
    *VC_GEE_pxq_divby_px1(MATRIX **, MATRIX **),
    /*multiple every element of matrix by the scalar, destroy the parameter matrix as needed, return a temporal new matrix*/
    *VC_GEE_scalar_times_matrix(double, MATRIX **),
    /*identity matrix of m*m */
    *VC_GEE_ident(int),
    /*column vector to diagonal vector*/
    *VC_GEE_form_diag(MATRIX **),
    /*subset of a matrix up to m-th row and n-th col from the left top corner*/
    *VC_GEE_corner(MATRIX **, int, int),
    /*row stacked lagged convariace matrix */
    *VC_GEE_covlag(MATRIX **, int, int),
    /*create a symmetric toeplitz matrix from either column or row vector, A[i,j] = A[i-1, j-1] for any i,j>0*/
    *VC_GEE_toeplitz(MATRIX **),
    /*set the rest of matrix equal to 0 except the band from the diagonal*/
    *VC_GEE_band(MATRIX **, int),
    /*extract m-th column to n-th column from the matrix*/
    *VC_GEE_extract_cols(MATRIX **, int, int),
    /*elmentwise pdf of the matrix*/
    *VC_GEE_matnpdf(MATRIX **),
    /*elmentwise cdf of the matrix*/
    *VC_GEE_matncdf(MATRIX **),
    /*elmentwise 1-exp(-exp(x)) of the matrix*/
    *VC_GEE_matanticlog(MATRIX **);
MATRIX *VC_GEE_diag_as_vec(MATRIX **);
MATRIX *VC_GEE_matsqrt(MATRIX **);
MATRIX *VC_GEE_mat1over(MATRIX **);
/*lump max and sum of the matrix*/
double VC_GEE_matmax(MATRIX **), VC_GEE_elsum(MATRIX **);

/*print, plug, free*/
int VC_GEE_matdump(MATRIX **),
    VC_GEE_plug(MATRIX **, MATRIX **, int, int);

int VC_GEE_split(MATRIX **, MATRIX **, MATRIX **[]),
    VC_GEE_nchanges(MATRIX **);

int Cgee(
    MATRIX **xin,
    MATRIX **yin,
    MATRIX **idin,
    MATRIX **nin,
    MATRIX **offsetin,
    int *nobs,
    int *p,
    int *parmvec,
    int *M_parm,
    MATRIX **betain,  // keep
    MATRIX **naivvar, // keep
    MATRIX **robvar,  // keep
    MATRIX **Rin,     // keep
    double *S_phi,   // keep
    double *tol,
    int *maxsz,  // keep
    int *S_iter, // keep
    int *silent,
    int *scale_fix,
    int *compatflag);

#ifndef MAX_NUM_CLUSTS
#define MAX_NUM_CLUSTS 5000
#endif

#define TOO_MANY_CLUSTS 100
#define UNKNOWN_LINK 101
#define UNKNOWN_VAR_MEAN_REL 102
#define M_DEP_SIZE_LT_M 103
#define MAXITER_EXCEEDED 104
#define LOGISTIC_DIVERGENCE 105

enum link_type
{
    VC_GEE_identity,
    logarithm,
    logit,
    reciprocal,
    probit,
    cloglog
};

enum var_mean_rel_type
{
    Gaussian,
    Poisson,
    Binomial,
    Gamma
};

enum corstruct_type
{
    independence,
    fixed,
    stat_M_dep,
    non_stat_M_dep,
    exchangeable,
    AR_M,
    unstructured
};

// Function to convert MATRIX to gsl_matrix
gsl_matrix *matrix_to_gsl(MATRIX **mat);

// Function to convert gsl_matrix to MATRIX
MATRIX *gsl_to_matrix(gsl_matrix *gslMat);

// Function to find the inverse of a MATRIX using GSL
MATRIX *matrix_inverse(MATRIX **mat);

// Print matrix out
void PrintMatrix(MATRIX **mat);

#endif