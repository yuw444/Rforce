#include "pvalue.h"

// Function to calculate the adjusted p-values
double *PAdjust(double *pvalues, size_t n, char method[])
{

  double *out = calloc(n, sizeof(double));
  if (out == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }

  if (n <= 1)
  {
    memcpy(out, pvalues, n * sizeof(double));
    return out;
  }

  if (n == 2 && strcmp(method, "hommel") == 0)
  {
    method = "hochberg";
  }

  if (strcmp(method, "bonferroni") == 0)
  {
    for (size_t i = 0; i < n; i++)
    {
      out[i] = pvalues[i] * n < 1 ? pvalues[i] * n : 1;
    }
    return out;
  }

  if (strcmp(method, "holm") == 0)
  {
    int *i = GetSeqInt(1L, n, 1L);
    int *o = Order(pvalues, n, false, true);
    int *ro = Order(o, n, false, false);
    double *pvaluesOrdered = calloc(n, sizeof(double));
    double * ptemp = calloc(n, sizeof(double));
    if (pvaluesOrdered == NULL || ptemp == NULL)
    {
      PRINT_LOCATION();
      printf("Memory allocation failed\n");
      exit(1);
    }
    for (size_t j = 0; j < n; j++)
    {
      pvaluesOrdered[j] = pvalues[o[j]];
    }

    double *n1_i_p = (double *)calloc(n, sizeof(double));
    if (n1_i_p == NULL)
    {
      PRINT_LOCATION();
      printf("Memory allocation failed\n");
      exit(1);
    }

    for (size_t j = 0; j < n; j++)
    {
      n1_i_p[j] = (n - i[j] + 1L) * pvaluesOrdered[j];
    }
    double *n1_i_p_cummax = Cummax(n1_i_p, n);
    for (size_t j = 0; j < n; j++)
    {
      ptemp[j] = n1_i_p_cummax[j] < 1 ? n1_i_p_cummax[j] : 1;
    }

    for (size_t j = 0; j < n; j++)
    {
      out[j] = ptemp[ro[j]];
    }

    Free(i);
    Free(o);
    Free(ro);
    Free(pvaluesOrdered);
    Free(n1_i_p);
    Free(n1_i_p_cummax);
    Free(ptemp);
    return out;
  }

  if (strcmp(method, "hommel") == 0)
  {

    int *i = GetSeqInt(1L, n, 1L);
    int *o = Order(pvalues, n, false, true);
    int *ro = Order(o, n, false, false);
    double *pvaluesOrdered = calloc(n, sizeof(double));
    double *ptemp = calloc(n, sizeof(double));
    double *npi = (double *)calloc(n, sizeof(double));
    double *q = (double *)calloc(n, sizeof(double));
    double *pa = (double *)calloc(n, sizeof(double));
    if (pvaluesOrdered == NULL || ptemp == NULL || npi == NULL || q == NULL || pa == NULL)
    {
      PRINT_LOCATION();
      printf("Memory allocation failed\n");
      exit(1);
    }
    for (size_t j = 0; j < n; j++)
    {
      pvaluesOrdered[j] = pvalues[o[j]];
    }

    for (size_t j = 0; j < n; j++)
    {
      npi[j] = n * pvaluesOrdered[j] / i[j];
    }
    double npi_min = Min(npi, n);
    for (size_t j = 0; j < n; j++)
    {
      q[j] = npi_min;
      pa[j] = npi_min;
    }
    for (size_t j = n - 1L; j >= 2; j--)
    {
      int *ij = GetSeqInt(1L, n - j + 1L, 1L);
      int *i2 = GetSeqInt(n - j + 2L, n, 1L);
      int *j2 = GetSeqInt(2L, j, 1L);
      double *jpi22j = (double *)calloc(j - 1L, sizeof(double));
      if (jpi22j == NULL)
      {
        PRINT_LOCATION();
        printf("Memory allocation failed\n");
        exit(1);
      }
      for (size_t k = 0; k < (j - 1L); k++)
      {
        jpi22j[k] = j * pvaluesOrdered[i2[k] - 1L] / j2[k];
      }

      double q1 = Min(jpi22j, j - 1L);
      
      for (size_t k = 0; k < n - j + 1L; k++)
      {
        q[ij[k] - 1L] = (j * pvaluesOrdered[ij[k] - 1L]) < q1 ? (j * pvaluesOrdered[ij[k] - 1]) : q1;
      }
      for (size_t k = 0; k < (j - 1L); k++)
      {
        q[i2[k] - 1L] = q[n - j];
      }
      for (size_t k = 0; k < n; k++)
      {
        pa[k] = pa[k] < q[k] ? q[k] : pa[k];
      }
      Free(ij);
      Free(i2);
      Free(j2);
      Free(jpi22j);
    }

    for (size_t j = 0; j < n; j++)
    {
      ptemp[j] = pa[j] > pvaluesOrdered[j] ? pa[j] : pvaluesOrdered[j];
    }

    for (size_t j = 0; j < n; j++)
    {
      out[j] = ptemp[ro[j]];
    }
    Free(i);
    Free(o);
    Free(ro);
    Free(pvaluesOrdered);
    Free(npi);
    Free(q);
    Free(pa);
    Free(ptemp);
    return out;
  }

  int *ir = GetSeqInt((int)n, 1L, -1L);
  int *or = Order(pvalues, n, true, true);
  int *ror = Order(or, n, false, false);
  double *pvaluesOrdered = calloc(n, sizeof(double));
  double *ptemp = calloc(n, sizeof(double));
  if (pvaluesOrdered == NULL || ptemp == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }
  for (size_t j = 0; j < n; j++)
  {
    pvaluesOrdered[j] = pvalues[or[j]];
  }

  if (strcmp(method, "hochberg") == 0)
  {
    for (size_t j = 0; j < n; j++)
    {
      ptemp[j] = (n + 1L - ir[j]) * pvaluesOrdered[j] < 1 ? (n + 1L - ir[j]) * pvaluesOrdered[j] : 1;
    }

    for (size_t j = 0; j < n; j++)
    {
      out[j] = ptemp[ror[j]];
    }
    Free(ir);
    Free(or);
    Free(ror);
    Free(pvaluesOrdered);
    Free(ptemp);
    return out;
  }

  if (strcmp(method, "BH") == 0)
  {
    for (size_t j = 0; j < n; j++)
    {
      ptemp[j] = n * 1.f / ir[j] * pvaluesOrdered[j];
    }

    double *ptemp1 = Cummin(ptemp, n);

    for (size_t j = 0; j < n; j++)
    {
      ptemp1[j] = 1.0 < ptemp1[j] ? 1.0 : ptemp1[j];
    }

    for (size_t j = 0; j < n; j++)
    {
      out[j] = ptemp1[ror[j]];
    }
    Free(ir);
    Free(or);
    Free(ror);
    Free(pvaluesOrdered);
    Free(ptemp);
    Free(ptemp1);
    return out;
  }

  if (strcmp(method, "BY") == 0)
  {
    double qsum = 0.f;
    for (size_t j = 0; j < n; j++)
    {
      qsum += 1.f / ir[j];
    }
    for (size_t j = 0; j < n; j++)
    {
      ptemp[j] = qsum * (n * 1.f) / ir[j] * pvaluesOrdered[j];
    }
    double *ptemp1 = Cummin(ptemp, n);

    for (size_t j = 0; j < n; j++)
    {
      ptemp1[j] = 1.0 < ptemp1[j] ? 1.0 : ptemp1[j];
    }

    for (size_t j = 0; j < n; j++)
    {
      out[j] = ptemp1[ror[j]];
    }
    Free(ir);
    Free(or);
    Free(ror);
    Free(pvaluesOrdered);
    Free(ptemp);
    Free(ptemp1);
    return out;
  }

  Free(out);
  return NULL;
}