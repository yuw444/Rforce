/*
@author : Yu Wang
*/

#ifndef UTILS_H
#define UTILS_H
#define NA 999999.0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stddef.h>
#include <stdbool.h>
#include <omp.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#define BUFFER 1024 * 1024

#define PRINT_LOCATION() printf("\nFile: %s:%d : ", __FILE__, __LINE__ - 3)
// #define FREE(p) if(*p != NULL) {free(*p); *p = NULL;}

// Macro to check for NA
#define IS_NA(x) ((x) == NA)

// struct for double elements and its length
typedef struct
{
  size_t nElements;
  double *Elements;
} ElementStruct;

// print an array
void PrintArrayInt(int *array, size_t n);
void PrintArrayDouble(double *array, size_t n);

double **Allocate2DArray (size_t nrow, size_t ncol);
double **Copy2DArray (double **array, size_t nrow, size_t ncol);
void Free2DArray (double **array, size_t nrow);

int **Allocate2DArrayInt(size_t nrow, size_t ncol);
void Free2DArrayInt(int **array, size_t nrow);

double ***Allocate3DArray (size_t nlayers, size_t nrow, size_t ncol);
void Free3DArray (double ***array, size_t nlayers, size_t nrow);

int *GetSeqInt(int start, int end, int step);

// get a double sequence
double *GetSeqDouble(double start, double end, double step);

// sample k elements with or without replacement from a integer array of length n using seed
int *SampleInt(int *arrayIn, size_t n, size_t nSample, unsigned int replace, unsigned int seed);

// NA remover for a double array with n elements
ElementStruct *RemoveNA(double *arrayIn, size_t n, double NA_value);

// sample k elements with or without replacement from a double array of length n using seed
double *SampleDouble(double *arrayIn, size_t n, size_t nSample, unsigned int replace, unsigned int seed);

// used for qsort
int vsD(const void *a, const void *b);

// used for qsort
int vsI(const void *a, const void *b);

// get unique elements and length of a double array with n elements
ElementStruct *Unique(double *arrayIn, size_t n);

// get transpose of a double array with nrows and ncols
double **Transpose(double **arrayIn, size_t nrows, size_t ncols);

// get nth column of a double array with nrows and ncols
double *GetCol(double **arrayIn, size_t nrows, size_t ncols, unsigned int nthCol);

// get quantiles of a double array with n elements
double *Quantile(double *arrayIn, size_t n, double *quantiles, size_t nQuantiles);

// get split candidates for a double array with n elements with nsplits
// a wrapper for Unique
ElementStruct *SplitCandidates(double *arrayIn, size_t n, size_t nsplits);

// get number of row of a csv file
size_t GetNrowCSV(char *filename);

// get number of column of a csv file
size_t GetNcolCSV(char *filename);

// read a csv file into a 2D double array with nrows and ncols using fscanf
void ReadCSV(char *filename, double ***data, unsigned int header);
void ReadCSVInt(char *filename, int ***data, unsigned int header);


// write a 2D double/int array with nrows and ncols into a csv file
void WriteCSV(double **data, FILE *out, size_t nrows, size_t ncols);
void WriteCSVInt(int **data, FILE *out, size_t nrows, size_t ncols);

// get sum of double array
double Sum(double *arrayIn, size_t n);

// get max of double array
double Max(double *arrayIn, size_t n);
double Min(double *arrayIn, size_t n);

// get mean and row mean of double array;
double Mean(double *arrayIn, size_t n);
double *RowMean(double **arrayIn, size_t nrows, size_t ncols);
double *RowMeanExcludingZeros(double **arrayIn, size_t nrows, size_t ncols); // for OOB prediction
double NthColMeanExcluding(double **arrayIn, size_t nrows, unsigned int nthCol, double exclude); // for dynamic imputation

// get variance of double array
double Var(double *arrayIn, size_t n);

// get median of double array
double Median(double *arrayIn, size_t n);

// delete non-empty directory
int DeleteDir (const char *path);

// get table of unique elements and their counts from an integer array
void TableInt(int *arrayIn, size_t n, int ***tableOut, size_t *nUniques);

typedef struct {
  double value;
  int index;
} IndexedValue;

// compare IndexedValue by value for qsort
// increasing
int vsP(const void *a, const void *b);
// decreasing
int vsPI(const void *a, const void *b);

int *Order(void *arrayIn, size_t n, bool decreasing, bool is_double);
double *Cummax(double *arrayIn, size_t n);
double *Cummin(double *arrayIn, size_t n);

double *Permute(double *arrayIn, size_t n, unsigned int seed);
double **ColsPermute(double **arrayIn, size_t nrows, size_t ncols, unsigned int *colsToPermute, size_t ncolsToPermute, unsigned int seed);

void MkdirRecursive(const char *path);
#endif