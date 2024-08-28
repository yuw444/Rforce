#include "split.h"
#include "surv.h"

// note that the last two columns of auxiliaryFeatures are X and status
double *LeafOutputIntervalDynamic(
    size_t nrow,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIU,
    size_t nUnits,
    size_t lenOutput)
{

  // get X and status from auxiliaryFeatures
  double *X = GetCol(auxiliaryFeatures, nrow, ncolsAuxiliary, ncolsAuxiliary-2);
  double *status = GetCol(auxiliaryFeatures, nrow, ncolsAuxiliary, ncolsAuxiliary-1);

  double *statusInverse = (double *)malloc(nrow * sizeof(double));

  for (int i = 0; i < nrow; i++)
  {
    if (fabs(status[i]) < 1e-9)
    {
      statusInverse[i] = 1.0;
    }
    else
    {
      statusInverse[i] = 0.0;
    }
  }

  // calculate the breakpoint of intervals
  double *breakpoints = (double *)malloc((nUnits + 1) * sizeof(double));
  breakpoints[0] = unitsOfCPIU[0];

  for (int i = 0; i < nUnits; i++)
  {
    breakpoints[i+1] = breakpoints[i] + unitsOfCPIU[i];
  }

  // fit kaplan-meier curve w.r.t censoring
  KMResult *Gt = km(nrow, X, statusInverse);
  KMResult *wtPatient; // placeholder for inverse probability of censoring weight w.r.t time of interest

  // get pseudo risktime for each patient each interval
  double **pseudoRiskTimes = Allocate2DArray(nrow, nUnits);

  for (int i = 0; i < nrow; i++)
  {
    wtPatient = wt(Gt, X[i], statusInverse[i]);
    for (int j = 0; j < nUnits; j++)
    {
      pseudoRiskTimes[i][j] = pseudoRiskTime(wtPatient, breakpoints[j], breakpoints[j+1]);
    }
    freeKMResult(wtPatient);
  }

  // get the output
  double *output = (double *)calloc(lenOutput, sizeof(double));
  double *Y = (double *)calloc(nrow, sizeof(double));
  double *rt = (double *)calloc(nrow, sizeof(double));

  for (int i = 0; i < nrow; i++)
  {
    for (int j = 0; j < nUnits; j++)
    {
      Y[j] += designMatrixY[i][ncolsDesign + j];
      rt[j] += pseudoRiskTimes[i][j];
    }
  }

  for (int i = 0; i < nUnits; i++)
  {
    output[i] = Y[i] / (rt[i] + 1e-9);
  } 

  free(Y);
  free(rt);
  free(X);
  free(status);
  free(statusInverse);
  free(breakpoints);
  freeKMResult(Gt);
  free2DArray(pseudoRiskTimes, nrow);

  return output;

}