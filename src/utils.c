#include "utils.h"

const double NA_DOUBLE = 0.0 / 0.0; // define NA as NaN

void PrintArrayInt(int *array, size_t n)
{
  if (array == NULL)
  {
    printf("NULL\n");
    return;
  }
  for (int i = 0; i < n; i++)
  {
    printf("%d ", array[i]);
  }
  printf("\n");
}

void PrintArrayDouble(double *array, size_t n)
{
  if (array == NULL)
  {
    printf("NULL\n");
    return;
  }
  for (int i = 0; i < n; i++)
  {
    if (isnan(array[i]))
      printf("NA ");
    else
      printf("%.3f ", array[i]);
  }
  printf("\n");
}

double **Allocate2DArray(size_t nrow, size_t ncol)
{
  double **array = (double **)malloc(nrow * sizeof(double *));
  if (array == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  for (int i = 0; i < nrow; i++)
  {
    array[i] = (double *)calloc(ncol, sizeof(double));
    if (array[i] == NULL)
    {
      PRINT_LOCATION();
      printf("Memory allocation failed\n");
      exit(1);
    }
  }
  return array;
}

double **Copy2DArray(double **array, size_t nrow, size_t ncol)
{
  double **arrayCopy = Allocate2DArray(nrow, ncol);
  for (int i = 0; i < nrow; i++)
  {
    for (int j = 0; j < ncol; j++)
    {
      arrayCopy[i][j] = array[i][j];
    }
  }
  return arrayCopy;
}

int **Allocate2DArrayInt(size_t nrow, size_t ncol)
{
  int **array = (int **)malloc(nrow * sizeof(int *));
  if (array == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  for (int i = 0; i < nrow; i++)
  {
    array[i] = (int *)calloc(ncol, sizeof(int));
    if (array[i] == NULL)
    {
      PRINT_LOCATION();
      printf("Memory allocation failed\n");
      exit(1);
    }
  }
  return array;
}

void Free2DArray(double **array, size_t nrow)
{
  if (array == NULL)
    return;
  for (int i = 0; i < nrow; i++)
  {

    free(array[i]);
  }
  free(array);
}

void Free2DArrayInt(int **array, size_t nrow)
{
  if (array == NULL)
    return;
  for (int i = 0; i < nrow; i++)
  {

    free(array[i]);
  }
  free(array);
}

double ***Allocate3DArray(size_t nlayers, size_t nrow, size_t ncol)
{
  double ***array = (double ***)malloc(nlayers * sizeof(double **));
  if (array == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  for (int i = 0; i < nlayers; i++)
  {
    array[i] = Allocate2DArray(nrow, ncol);
    if (array[i] == NULL)
    {
      PRINT_LOCATION();
      printf("Memory allocation failed\n");
      exit(1);
    }
  }
  return array;
}

void Free3DArray(double ***array, size_t nlayers, size_t nrow)
{
  if (array == NULL)
    return;
  for (int i = 0; i < nlayers; i++)
  {
    Free2DArray(array[i], nrow);
  }
  free(array);
}

int *GetSeqInt(int start, int end, int step)
{
  if (step == 0)
  {
    PRINT_LOCATION();
    printf("Error: step cannot be 0 in `GetSeqInt()`\n");
    exit(1);
  }

  if ((end - start) * step < 0)
  {
    PRINT_LOCATION();
    printf("Error: step must be positive when start < end, or negative when start > end in `GetSeqInt()`\n");
    exit(1);
  }
  size_t n = (end - start) / step + 1;
  int *seq = (int *)calloc(n, sizeof(int));
  if (seq == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  for (int i = 0; i < n; i++)
  {
    seq[i] = start + i * step;
  }
  return seq;
}

double *GetSeqDouble(double start, double end, double step)
{
  if (step == 0.0)
  {
    PRINT_LOCATION();
    printf("Error: step cannot be 0 in `GetSeqDouble()`\n");
    exit(1);
  }

  if ((end - start) * step < 0.0)
  {
    PRINT_LOCATION();
    printf("Error: step must be positive when start < end, or negative when start > end in `GetSeqDouble()`\n");
    exit(1);
  }

  int n = (int)((end - start) / step + 1);
  double *seq = (double *)calloc(n, sizeof(double));
  if (seq == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  for (int i = 0; i < n; i++)
  {
    seq[i] = start + i * step;
  }
  return seq;
}

int *SampleInt(int *arrayIn, size_t n, size_t nSample, unsigned int replace, unsigned int seed)
{
  int *sampleOut = calloc(nSample, sizeof(int));
  if (sampleOut == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }

  /* local RNG state seeded from the seed argument */
  xoshiro256pp_state rng;
  xoshiro256pp_seed(&rng, (uint64_t)seed);

  if (replace == 0)
  {
    if (nSample > n)
    {
      free(sampleOut);
      PRINT_LOCATION();
      printf("Sample size must be smaller than population size when sampling without replacement in `SampleInt()`.\n");
      exit(1);
    }

    int *arrayInCopy = calloc(n, sizeof(int));
    if (arrayInCopy == NULL)
    {
      free(sampleOut);
      PRINT_LOCATION();
      printf("Memory allocation failed\n");
      exit(1);
    }
    memcpy(arrayInCopy, arrayIn, n * sizeof(int));

    size_t m = n;
    for (size_t i = 0; i < nSample; i++)
    {
      size_t index = xoshiro256pp_uniform_index(&rng, m);
      sampleOut[i] = arrayInCopy[index];
      if (index != (m - 1))
      {
        arrayInCopy[index] = arrayInCopy[m - 1];
      }
      m--;
    }

    free(arrayInCopy);
  }
  else
  {
    for (size_t i = 0; i < nSample; i++)
    {
      size_t index = xoshiro256pp_uniform_index(&rng, n);
      sampleOut[i] = arrayIn[index];
    }
  }

  return sampleOut;
}

double *SampleDouble(double *arrayIn, size_t n, size_t nSample, unsigned int replace, unsigned int seed)
{
  double *sampleOut = calloc(nSample, sizeof(double));
  if (sampleOut == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }

  xoshiro256pp_state rng;
  xoshiro256pp_seed(&rng, (uint64_t)seed);

  if (replace == 0)
  {
    if (nSample > n)
    {
      free(sampleOut);
      PRINT_LOCATION();
      printf("Sample size must be smaller than population size when sampling without replacement in `SampleDouble()`.\n");
      exit(1);
    }

    double *arrayInCopy = calloc(n, sizeof(double));
    if (arrayInCopy == NULL)
    {
      free(sampleOut);
      PRINT_LOCATION();
      printf("Memory allocation failed\n");
      exit(1);
    }
    memcpy(arrayInCopy, arrayIn, n * sizeof(double));

    size_t m = n;
    for (size_t i = 0; i < nSample; i++)
    {
      size_t index = xoshiro256pp_uniform_index(&rng, m);
      sampleOut[i] = arrayInCopy[index];
      if (index != (m - 1))
      {
        arrayInCopy[index] = arrayInCopy[m - 1];
      }
      m--;
    }

    free(arrayInCopy);
  }
  else
  {
    for (size_t i = 0; i < nSample; i++)
    {
      size_t index = xoshiro256pp_uniform_index(&rng, n);
      sampleOut[i] = arrayIn[index];
    }
  }

  return sampleOut;
}

// NA remover for a double array with n elements
ElementStruct *RemoveNA(double *arrayIn, size_t n)
{
  ElementStruct *out = malloc(sizeof(ElementStruct));
  if (out == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  double *arrayOut = calloc(n, sizeof(double));
  if (arrayOut == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  size_t j = 0;
  for (size_t i = 0; i < n; i++)
  {
    if (!isnan(arrayIn[i]))
    {
      arrayOut[j++] = arrayIn[i];
    }
  }
  arrayOut = (double *)realloc(arrayOut, j * sizeof(double));
  *out = (ElementStruct){j, arrayOut};
  return out;
}

// used for qsort
int vsD(const void *a, const void *b)
{
  if (*(double *)a > *(double *)b)
    return 1;
  else if (*(double *)a < *(double *)b)
    return -1;
  else
    return 0;
}

// used for qsort
int vsI(const void *a, const void *b)
{
  return (*(int *)a - *(int *)b);
}

ElementStruct *Unique(double *arrayIn, size_t n)
{
  // in case NA is present, safely remove them first
  ElementStruct *naRemoved = RemoveNA(arrayIn, n);
  // PrintArrayDouble(naRemoved->Elements, naRemoved->nElements);

  qsort(naRemoved->Elements, naRemoved->nElements, sizeof(double), vsD);

  double *uniqueElements = (double *)calloc(naRemoved->nElements, sizeof(double));
  if (uniqueElements == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  size_t nUniques = 1;
  uniqueElements[0] = naRemoved->Elements[0];

  for (int i = 1; i < naRemoved->nElements; ++i)
  {
    // make sure it is not a missing value
    if (fabs((naRemoved->Elements)[i - 1] - (naRemoved->Elements)[i]) > 1e-9)
    {
      uniqueElements[nUniques++] = naRemoved->Elements[i];
    }
  }

  uniqueElements = (double *)realloc(uniqueElements, nUniques * sizeof(double));

  ElementStruct *out = malloc(sizeof(ElementStruct));
  if (out == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }

  *out = (ElementStruct){nUniques, uniqueElements};

  free(naRemoved->Elements);
  free(naRemoved);

  return out;
}

double **Transpose(double **arrayIn, size_t nrows, size_t ncols)
{
  double **arrayOut = Allocate2DArray(ncols, nrows);

  for (int i = 0; i < nrows; i++)
  {
    for (int j = 0; j < ncols; j++)
    {
      arrayOut[j][i] = arrayIn[i][j];
    }
  }

  return arrayOut;
}

double *GetCol(double **arrayIn, size_t nrows, size_t ncols, unsigned int nthCol)
{
  if (nthCol > ncols)
  {
    PRINT_LOCATION();
    printf("Column index out of bounds in `GetCol()`.\n");
    exit(1);
  }

  double *colOut = calloc(nrows, sizeof(double));
  if (colOut == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }

  for (int i = 0; i < nrows; i++)
  {
    colOut[i] = arrayIn[i][nthCol];
  }

  return colOut;
}

double *Quantile(double *arrayIn, size_t n, double *quantiles, size_t nQuantiles)
{
  double *arrayInCopy = calloc(n, sizeof(double));
  double *quantileOut = calloc(nQuantiles, sizeof(double));
  if (arrayInCopy == NULL || quantileOut == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }

  memcpy(arrayInCopy, arrayIn, n * sizeof(double));

  qsort(arrayInCopy, n, sizeof(double), vsD);

  for (int i = 0; i < nQuantiles; i++)
  {
    quantileOut[i] = arrayInCopy[(int)(quantiles[i] * (n - 1))];
  }

  free(arrayInCopy);

  return quantileOut;
}

ElementStruct *SplitCandidates(double *arrayIn, size_t n, size_t nsplits)
{
  ElementStruct *uniques = Unique(arrayIn, n);
  // PrintArrayDouble(uniques->Elements, uniques->nElements);

  if (nsplits >= uniques->nElements)
  {
    return uniques;
  }
  else
  {

    // free(uniques->Elements);

    // free(uniques);

    ElementStruct *quantileOut = malloc(sizeof(ElementStruct));
    if (quantileOut == NULL)
    {
      PRINT_LOCATION();
      printf("Memory allocation failed\n");
      exit(1);
    }
    // old method that use quantiles as the split candidates
    // double *quantiles = GetSeqDouble(1.0 / (1 + nsplits), 1.0 - 1e-6, 1.0 / (1 + nsplits));
    // printf("nsplits: %ld; nuniques: %ld\n", nsplits, uniques->nElements);
    *quantileOut = (ElementStruct){nsplits, SampleDouble(uniques->Elements, uniques->nElements, nsplits, 0L, 926L)};
    free(uniques->Elements);
    free(uniques);
    // free(quantiles);
    return quantileOut;
  }
}

size_t GetNrowCSV(char *filename)
{
  char line[BUFFER];
  int nrows = 0;
  FILE *stream = fopen(filename, "r");
  if (stream == NULL)
  {
    PRINT_LOCATION();
    printf("Error: cannot open file %s in `GetNrowCSV`.\n", filename);
    exit(1);
  }
  while (fgets(line, BUFFER, stream))
  {
    nrows++;
  }
  fclose(stream);
  return nrows;
}

size_t GetNcolCSV(char *filename)
{
  char line[BUFFER];
  int ncols = 0;
  FILE *stream = fopen(filename, "r");
  if (stream == NULL)
  {
    PRINT_LOCATION();
    printf("Error: cannot open file %s\n", filename);
    exit(1);
  }
  fgets(line, BUFFER, stream);
  char *token = strtok(line, ",");
  while (token != NULL)
  {
    ncols++;
    token = strtok(NULL, ",");
  }
  fclose(stream);
  return ncols;
}

// read csv file to data 2d array
void ReadCSV(char *filename, double ***data, unsigned int header)
{
  char line[BUFFER];

  FILE *stream = fopen(filename, "r");

  if (stream == NULL)
  {
    PRINT_LOCATION();
    printf("Error: open %s in `ReadCSV()`.\n", filename);
    exit(1);
  }
  if (header == 1)
  {
    fgets(line, BUFFER, stream);
  }

  int i = 0;
  while (fgets(line, BUFFER, stream))
  {
    int j = 0;
    char *tok;
    char *tmp = strdup(line);
    for (tok = strtok(line, ",");
         tok && *tok;
         j++, tok = strtok(NULL, ",\n"))
    {
      if (strcmp(tok, "NA") == 0 || strcmp(tok, ".") == 0)
      {
        (*data)[i][j] = NA_DOUBLE;
      }
      else
      {
        (*data)[i][j] = atof(tok);
      }
    }
    i++;

    free(tmp);
  }
  fclose(stream);
}

// read csv file to data 2d array
void ReadCSVInt(char *filename, int ***data, unsigned int header)
{
  char line[BUFFER];

  FILE *stream = fopen(filename, "r");
  if (stream == NULL)
  {
    PRINT_LOCATION();
    printf("Error: open %s\n", filename);
    exit(1);
  }
  if (header == 1)
  {
    fgets(line, BUFFER, stream);
  }

  int i = 0;
  while (fgets(line, BUFFER, stream))
  {
    int j = 0;
    char *tok;
    char *tmp = strdup(line);
    for (tok = strtok(line, ",");
         tok && *tok;
         j++, tok = strtok(NULL, ",\n"))
    {
      (*data)[i][j] = atoi(tok);
      // printf("i: %d, j: %d\n", i, j);
    }
    i++;

    free(tmp);
  }
  fclose(stream);
}

void WriteCSV(double **data, FILE *out, size_t nrows, size_t ncols)
{
  for (int i = 0; i < nrows; i++)
  {
    for (int j = 0; j < ncols - 1; j++)
    {
      fprintf(out, "%.5f,", data[i][j]);
    }
    fprintf(out, "%.5f\n", data[i][ncols - 1]);
  }
}

void WriteCSVInt(int **data, FILE *out, size_t nrows, size_t ncols)
{
  for (int i = 0; i < nrows; i++)
  {
    for (int j = 0; j < ncols - 1; j++)
    {
      fprintf(out, "%d,", data[i][j]);
    }
    fprintf(out, "%d\n", data[i][ncols - 1]);
  }
}

double Sum(double *arrayIn, size_t n)
{
  double sum = 0.0;
  for (int i = 0; i < n; i++)
  {
    sum += arrayIn[i];
  }
  return sum;
}

double Mean(double *arrayIn, size_t n)
{
  return Sum(arrayIn, n) / n;
}

double *RowMean(double **arrayIn, size_t nrows, size_t ncols)
{
  double *rowMean = calloc(nrows, sizeof(double));
  if (rowMean == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  for (int i = 0; i < nrows; i++)
  {
    rowMean[i] = Mean(arrayIn[i], ncols);
  }
  return rowMean;
}

double Max(double *arrayIn, size_t n)
{
  double max = arrayIn[0];
  for (int i = 1; i < n; i++)
  {
    if (arrayIn[i] > max)
    {
      max = arrayIn[i];
    }
  }
  return max;
}

double Min(double *arrayIn, size_t n)
{
  double min = arrayIn[0];
  for (int i = 1; i < n; i++)
  {
    if (arrayIn[i] < min)
    {
      min = arrayIn[i];
    }
  }
  return min;
}

double *RowMeanExcludingZeros(double **arrayIn, size_t nrows, size_t ncols)
{
  double *rowMean = calloc(nrows, sizeof(double));
  if (rowMean == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  for (int i = 0; i < nrows; i++)
  {
    int n = 0;
    for (int j = 0; j < ncols; j++)
    {
      if (fabs(arrayIn[i][j] - 0) > 1e-9)
      {
        rowMean[i] += arrayIn[i][j];
        n++;
      }
    }
    rowMean[i] /= n;
  }
  return rowMean;
}

// the user should make sure that the nthCol is out of bounds
double NthColMean(
    double **arrayIn,
    size_t nrows,
    unsigned int nthCol)
{
  double tempSum = 0.0;
  double tempCounts = 1e-9;

  for (int i = 0; i < nrows; i++)
  {
    if (!isnan(arrayIn[i][nthCol]))
    {
      tempSum += arrayIn[i][nthCol];
      tempCounts++;
    }
  }

  return tempSum / tempCounts;
}

double Var(double *arrayIn, size_t n)
{
  double mean = Mean(arrayIn, n);
  double var = 0.0;
  for (int i = 0; i < n; i++)
  {
    var += pow(arrayIn[i] - mean, 2);
  }
  return var / (n - 1);
}

// get median of double array
double Median(double *arrayIn, size_t n)
{
  double *arrayInCopy = calloc(n, sizeof(double));
  if (arrayInCopy == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  memcpy(arrayInCopy, arrayIn, n * sizeof(double));

  qsort(arrayInCopy, n, sizeof(double), vsD);

  double median = 0.0;
  if (n % 2 == 0)
  {
    median = (arrayInCopy[n / 2 - 1] + arrayInCopy[n / 2]) / 2;
  }
  else
  {
    median = arrayInCopy[n / 2];
  }

  free(arrayInCopy);

  return median;
}

// get table of unique elements and their counts from an integer array
void TableInt(int *arrayIn, size_t n, int ***table, size_t *nUniques)
{
  int *arrayInCopy = calloc(n, sizeof(int));
  if (arrayInCopy == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  int **tableTemp = *table;

  memcpy(arrayInCopy, arrayIn, n * sizeof(int));

  qsort(arrayInCopy, n, sizeof(int), vsI);
  tableTemp[0][0] = arrayInCopy[0];
  tableTemp[1][0] = 1;

  *nUniques = 1;

  for (int i = 1; i < n; ++i)
  {
    if (arrayInCopy[i - 1] != arrayInCopy[i])
    {
      tableTemp[0][*nUniques] = arrayInCopy[i];
      tableTemp[1][*nUniques] = 1;
      (*nUniques)++;
    }
    else
    {
      tableTemp[1][(*nUniques) - 1]++;
    }
  }

  free(arrayInCopy);
  tableTemp[0] = (int *)realloc(tableTemp[0], (*nUniques) * sizeof(int));
  tableTemp[1] = (int *)realloc(tableTemp[1], (*nUniques) * sizeof(int));
}

int DeleteDir(const char *path)
{
  DIR *dir = opendir(path);
  struct dirent *entry;
  char filepath[1024];

  if (dir == NULL)
  {
    PRINT_LOCATION();
    perror("Error opening directory");
    return -1;
  }

  while ((entry = readdir(dir)) != NULL)
  {
    if (entry->d_type == DT_DIR)
    {
      if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
      {
        continue;
      }
      snprintf(filepath, sizeof(filepath), "%s/%s", path, entry->d_name);
      DeleteDir(filepath);
    }
    else
    {
      snprintf(filepath, sizeof(filepath), "%s/%s", path, entry->d_name);
      unlink(filepath);
    }
  }

  closedir(dir);

  if (rmdir(path) == -1)
  {
    PRINT_LOCATION();
    perror("Error deleting directory");
    return -1;
  }

  return 0;
}

int vsP(const void *a, const void *b)
{
  double x = ((IndexedValue *)a)->value;
  double y = ((IndexedValue *)b)->value;
  if (x < y)
    return -1;
  else
    return 1;
}

int vsPI(const void *a, const void *b)
{
  double x = ((IndexedValue *)a)->value;
  double y = ((IndexedValue *)b)->value;
  if (x < y)
    return 1;
  else
    return -1;
}

int *Order(void *arrayIn, size_t n, bool decreasing, bool is_double)
{
  IndexedValue *indexedValues = calloc(n, sizeof(IndexedValue));
  if (indexedValues == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }

  for (size_t i = 0; i < n; i++)
  {
    indexedValues[i].index = i;
    if (is_double)
    {
      indexedValues[i].value = ((double *)arrayIn)[i];
    }
    else
    {
      indexedValues[i].value = ((int *)arrayIn)[i];
    }
  }

  if (decreasing)
  {
    qsort(indexedValues, n, sizeof(IndexedValue), vsPI);
  }
  else
  {
    qsort(indexedValues, n, sizeof(IndexedValue), vsP);
  }

  int *order = calloc(n, sizeof(int));
  if (order == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  for (size_t i = 0; i < n; i++)
  {
    order[i] = indexedValues[i].index;
  }

  free(indexedValues);
  return order;
}

double *Cummax(double *arrayIn, size_t n)
{
  double *cummax = calloc(n, sizeof(double));
  if (cummax == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  cummax[0] = arrayIn[0];
  for (size_t i = 1; i < n; i++)
  {
    cummax[i] = arrayIn[i] > cummax[i - 1] ? arrayIn[i] : cummax[i - 1];
  }
  return cummax;
}

double *Cummin(double *arrayIn, size_t n)
{
  double *cummin = calloc(n, sizeof(double));
  if (cummin == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  cummin[0] = arrayIn[0];
  for (size_t i = 1; i < n; i++)
  {
    cummin[i] = arrayIn[i] < cummin[i - 1] ? arrayIn[i] : cummin[i - 1];
  }
  return cummin;
}

double *Permute(double *arrayIn, size_t n, unsigned int seed)
{
  double *arrayOut = calloc(n, sizeof(double));
  if (arrayOut == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  memcpy(arrayOut, arrayIn, n * sizeof(double));

  // Local RNG instead of srand/rand
  xoshiro256pp_state rng;
  xoshiro256pp_seed(&rng, (uint64_t)seed);

  for (size_t i = 0; i < n; i++)
  {
    size_t j = xoshiro256pp_uniform_index(&rng, n);
    double temp = arrayOut[i];
    arrayOut[i] = arrayOut[j];
    arrayOut[j] = temp;
  }

  return arrayOut;
}

double **ColsPermute(double **arrayIn, size_t nrows, size_t ncols, unsigned int *colsToPermute, size_t ncolsToPermute, unsigned int seed)
{
  double **arrayOut = Copy2DArray(arrayIn, nrows, ncols);

  for (size_t i = 0; i < ncolsToPermute; i++)
  {
    double *temp = GetCol(arrayIn, nrows, ncols, colsToPermute[i]);
    double *permutedCol = Permute(temp, nrows, seed + i);
    for (size_t j = 0; j < nrows; j++)
    {
      arrayOut[j][colsToPermute[i]] = permutedCol[j];
    }
    free(temp);
    free(permutedCol);
  }

  return arrayOut;
}

void MkdirRecursive(const char *path)
{
  char temp[256];
  char *pos = NULL;
  size_t len;

  // Copy the path to a temporary buffer
  strncpy(temp, path, sizeof(temp));
  len = strlen(temp);

  if (temp[len - 1] == '/')
  {
    temp[len - 1] = '\0';
  }

  // Iterate through the path and create directories
  for (pos = temp + 1; *pos; pos++)
  {
    if (*pos == '/')
    {
      *pos = '\0';
      if (mkdir(temp, 0700) != 0 && errno != EEXIST)
      {
        perror("mkdir");
        return;
      }
      *pos = '/';
    }
  }

  // Create the final directory
  if (mkdir(temp, 0700) != 0 && errno != EEXIST)
  {
    perror("mkdir");
  }
}