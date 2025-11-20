#include "split.h"

int _k = 4;

int _noPseudo = 0;
int _pseudoRisk1 = 0; // _pseudoRisk1 = 1; pseudo risk time estimate for the entire population
int _pseudoRisk2 = 0; // -pseudoRisk2 = 1; pseudo risk time estimate at each tree
int _dynamicRisk = 0; 

int _noPhi = 0; // _noPhi = 1; phi = 1 all the time
int _phi1 = 0;  // _phi1 = 1; phi is estimated for the entire population
int _phi2 = 0;  // _phi2 = 1; phi is estimated at each tree
int _dynamicPhi = 0;

int _longformat = 0;
int _gee = 0;
int _interaction = 0;
int _asympotic = 0;

int _verbose = 0;

double **_treePhi = NULL;
MATRIX *_Rin = NULL;
char *_padjust = NULL;

size_t *SplitSizes(
    unsigned int varIndex,      // var index
    double splitValue,          // best cutoff
    size_t nrows,               // nrows of design matrix, number of variables
    double **designMatrixY,     // design matrix (X) and response Y
    double **auxiliaryFeatures) // auxiliaryFeatures
{
  // split the dataset
  size_t *out = (size_t *)malloc(2 * sizeof(size_t));
  // imputed value for missing value
  double imputedMean = NthColMean(designMatrixY, nrows, varIndex);
  // spot for imputed value or actual value
  double tempValue = 0.0;

  for (size_t i = 0; i < nrows; i++)
  {
    // if the value is missing, impute it with the mean
    if (isnan(designMatrixY[i][varIndex]))
    {
      tempValue = imputedMean;
    }
    else
    {
      tempValue = designMatrixY[i][varIndex];
    }

    if (tempValue <= splitValue)
    {
      out[0]++;
    }
    else
    {
      out[2]++;
    }
  }

  return out;
}

DecisionTreeData *SplitDataset(
    unsigned int varIndex,      // var index
    double splitValue,          // best cutoff
    size_t nrows,               // nrows of design matrix, number of variables
    double **designMatrixY,     // design matrix (X) and response Y
    double **auxiliaryFeatures) // auxiliaryFeatures)
{
  double **designMatrixYLeft = (double **)malloc(nrows * sizeof(double *));
  double **designMatrixYRight = (double **)malloc(nrows * sizeof(double *));
  double **auxiliaryFeaturesLeft = (double **)malloc(nrows * sizeof(double *));
  double **auxiliaryFeaturesRight = (double **)malloc(nrows * sizeof(double *));

  size_t nLeft = 0;
  size_t nRight = 0;
  // imputed value for missing value
  double imputedMean = NthColMean(designMatrixY, nrows, varIndex);
  // spot for imputed value or actual value
  double tempValue = 0.0;

  for (size_t i = 0; i < nrows; i++)
  {
    // if the value is missing, impute it with the mean
    if (isnan(designMatrixY[i][varIndex]))
    {
      tempValue = imputedMean;
    }
    else
    {
      tempValue = designMatrixY[i][varIndex];
    }

    if (tempValue <= splitValue)
    {
      designMatrixYLeft[nLeft] = designMatrixY[i];
      auxiliaryFeaturesLeft[nLeft] = auxiliaryFeatures[i];
      nLeft++;
    }
    else
    {
      designMatrixYRight[nRight] = designMatrixY[i];
      auxiliaryFeaturesRight[nRight] = auxiliaryFeatures[i];
      nRight++;
    }
  }

  // realloc the memory for actual size
  designMatrixYLeft = (double **)realloc(designMatrixYLeft, nLeft * sizeof(double *));
  designMatrixYRight = (double **)realloc(designMatrixYRight, nRight * sizeof(double *));
  auxiliaryFeaturesLeft = (double **)realloc(auxiliaryFeaturesLeft, nLeft * sizeof(double *));
  auxiliaryFeaturesRight = (double **)realloc(auxiliaryFeaturesRight, nRight * sizeof(double *));

  // create a struct to store the data for two daughter nodes
  DecisionTreeData *dataSplits = (DecisionTreeData *)malloc(2 * sizeof(DecisionTreeData));

  dataSplits[0] = (DecisionTreeData){nLeft, designMatrixYLeft, auxiliaryFeaturesLeft};
  dataSplits[1] = (DecisionTreeData){nRight, designMatrixYRight, auxiliaryFeaturesRight};

  return dataSplits;
}

int **BagMatrix(
    double *patientIds,
    size_t nrowsDesign,
    size_t nTrees,
    unsigned int seed)
{
    srand(seed);
    // sample patient with replacement

    int **bagMatrix = Allocate2DArrayInt(nTrees, nrowsDesign);
    ElementStruct *uniquePatients = Unique(patientIds, nrowsDesign);

    for (int i = 0; i < nTrees; i++)
    {

        double *patientSample = SampleDouble(
            uniquePatients->Elements,
            uniquePatients->nElements,
            uniquePatients->nElements,
            1,
            seed + i);

        for (int j = 0; j < uniquePatients->nElements; j++)
        {
            for (int k = 0; k < nrowsDesign; k++)
            {
                if (fabs(patientSample[j] - patientIds[k]) < 1e-9)
                {
                    bagMatrix[i][k] += 1;
                }
            }
        }
        free(patientSample);
    }

    free(uniquePatients->Elements);
    free(uniquePatients);

    return bagMatrix;
}

DecisionTreeData *BootStrapSample(double **designMatrixY,
                                  double **auxiliaryFeatures,
                                  size_t nrowsDesign,
                                  int *bootVector // length nrowsDesign
)
{

    double **designMatrixYIn = malloc(nrowsDesign * sizeof(double *)); // bootstrap sample has up to nrowsDesign rows
    double **auxiliaryFeaturesIn = malloc(nrowsDesign * sizeof(double *));
    double **designMatrixYOOB = malloc(nrowsDesign * sizeof(double *));
    double **auxiliaryFeaturesOOB = malloc(nrowsDesign * sizeof(double *));

    int kIn = 0;
    int kOOB = 0;

    for (int i = 0; i < nrowsDesign; i++)
    {
        // the element of bootVector is at least 0, could be greater than 1
        if (bootVector[i] != 0)
        {
            for (int j = 0; j < bootVector[i]; j++)
            {
                designMatrixYIn[kIn] = designMatrixY[i];
                auxiliaryFeaturesIn[kIn] = auxiliaryFeatures[i];
                kIn++;
            }
        }
        else
        {
            designMatrixYOOB[kOOB] = designMatrixY[i];
            auxiliaryFeaturesOOB[kOOB] = auxiliaryFeatures[i];
            kOOB++;
        }
    }

    designMatrixYOOB = realloc(designMatrixYOOB, kOOB * sizeof(double *));
    auxiliaryFeaturesOOB = realloc(auxiliaryFeaturesOOB, kOOB * sizeof(double *));

    DecisionTreeData *bootStrapDataSet = malloc(2 * sizeof(DecisionTreeData));

    bootStrapDataSet[0] = (DecisionTreeData){kIn, designMatrixYIn, auxiliaryFeaturesIn};
    bootStrapDataSet[1] = (DecisionTreeData){kOOB, designMatrixYOOB, auxiliaryFeaturesOOB};

    return bootStrapDataSet;
}


double *LeafOutput(
    size_t nrow,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIU,
    long treeId,
    size_t nUnits,
    size_t lenOutput)
{
  double *output = (double *)malloc(lenOutput * sizeof(double));

  double Y = 0.0;  // number of composite events
  double rt = 0.0; // pseudo risk time

  for (size_t i = 0; i < nrow; i++)
  {
    Y += designMatrixY[i][ncolsDesign];
    rt += auxiliaryFeatures[i][ncolsAuxiliary - 1];
  }
  double alpha = 1.0 / (_k * _k);
  double beta = alpha / (Y / (rt + 1e-9));

  output[0] = (alpha + Y) / (beta + rt + 1e-9); // calculate the average rate
  return output;
}

// Interval Specific LeafOutput
double *LeafOutputInterval(
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIUs,
    long treeId,
    size_t nUnits,
    size_t lenOutput)
{
  double *output = (double *)malloc(lenOutput * sizeof(double));
  if (output == NULL)
  {
    PRINT_LOCATION();
    printf("Memory allocation failed\n");
    exit(1);
  }

  if (_dynamicRisk == 1)
  {
    /******************For dynamic pseudo-risk time estimation**************************/
    // pseudo risk time at current nodes
    // pseudo risk time at tree-level is handle in GrowTree() of tree.c
    // no pseudo risk time and population level pseudo risk time are handled by the user input
    double **pseudoRiskTime = Allocate2DArray(nrows, nUnits);
    /******************For dynamic pseudo-risk time estimation**************************/
    // get X and status from auxiliaryFeatures
    double *X = GetCol(auxiliaryFeatures, nrows, ncolsAuxiliary, ncolsAuxiliary - 2);
    double *status = GetCol(auxiliaryFeatures, nrows, ncolsAuxiliary, ncolsAuxiliary - 1);
    double *statusInverse = (double *)malloc(nrows * sizeof(double));

    if (statusInverse == NULL)
    {
      PRINT_LOCATION();
      printf("Memory allocation failed\n");
      exit(1);
    }

    for (int i = 0; i < nrows; i++)
    {
      if (fabs(status[i]) > 0.99)
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
    breakpoints[0] = 0;

    for (int i = 0; i < nUnits; i++)
    {
      breakpoints[i + 1] = breakpoints[i] + unitsOfCPIUs[i];
    }

    // fit kaplan-meier curve w.r.t censoring
    KMResult *Gt = KaplanMeier(X, statusInverse, nrows);

    free(statusInverse);

    // placeholder for inverse probability of censoring weight w.r.t time of interest
    KMResult *wtPatient;

    for (int i = 0; i < nrows; i++)
    {
      wtPatient = Wt(Gt, X[i], status[i]);
      for (int j = 0; j < nUnits; j++)
      {
        pseudoRiskTime[i][j] = PseudoRiskTime(wtPatient, breakpoints[j], breakpoints[j + 1]);
      }
      FreeKMResult(wtPatient);
    }
    FreeKMResult(Gt);

    free(X);

    free(status);

    double *Y = (double *)calloc(lenOutput, sizeof(double));  // number of composite events
    double *rt = (double *)calloc(lenOutput, sizeof(double)); // pseudo risk time
    double YTotal = 0.0;
    double rtTotal = 0.0;

    // calculate interval specific Y and rt
    for (size_t i = 0; i < nrows; i++)
    {
      for (size_t j = 0; j < lenOutput; j++)
      {
        Y[j] += designMatrixY[i][ncolsDesign + j];
        rt[j] += pseudoRiskTime[i][j];
        YTotal += designMatrixY[i][ncolsDesign + j];
        rtTotal += pseudoRiskTime[i][j];
      }
    }

    double alpha = 1.0 / (_k * _k);
    double beta = alpha / (YTotal / (rtTotal + 1e-9));

    // calculate the average rate
    for (size_t i = 0; i < lenOutput; i++)
    {
      output[i] = (alpha + Y[i]) / (beta + rt[i] + 1e-9);
    }

    free(Y);

    free(rt);

    free(breakpoints);
    Free2DArray(pseudoRiskTime, nrows);
    return output;
  }

  double *Y = (double *)calloc(lenOutput, sizeof(double));  // number of composite events
  double *rt = (double *)calloc(lenOutput, sizeof(double)); // pseudo risk time
  double YTotal = 0.0;
  double rtTotal = 0.0;

  // calculate interval specific Y and rt
  for (size_t i = 0; i < nrows; i++)
  {
    for (size_t j = 0; j < lenOutput; j++)
    {
      Y[j] += designMatrixY[i][ncolsDesign + j];
      rt[j] += auxiliaryFeatures[i][j + 1];
      YTotal += designMatrixY[i][ncolsDesign + j];
      rtTotal += auxiliaryFeatures[i][j + 1];
    }
  }

  double alpha = 1.0 / (_k * _k);
  double beta = alpha / (YTotal / (rtTotal + 1e-9));

  // calculate the average rate
  for (size_t i = 0; i < lenOutput; i++)
  {
    output[i] = (alpha + Y[i]) / (beta + rt[i] + 1e-9);
  }
  free(Y);
  free(rt);
  return output;
}

// the original RF-SLAM only
double *PoissonLikelihood(
    fpLeafOutput leafOutputFunction,
    unsigned int splitIndex,
    double splitValue,
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIUs,
    long treeId,
    size_t nUnits,
    size_t lenOutput,
    size_t minNodeSize)
{
  // split the dataset
  DecisionTreeData *dataSplits = SplitDataset(
      splitIndex,
      splitValue,
      nrows,
      designMatrixY,
      auxiliaryFeatures);

  // get the number of rows in left and right nodes
  size_t nLeft = dataSplits[0].nrows;
  size_t nRight = dataSplits[1].nrows;

  // check if the number of rows in left and right nodes are zero
  if (nLeft <= minNodeSize || nRight == minNodeSize)
  {
    free(dataSplits[0].designMatrixY);
    free(dataSplits[0].auxiliaryFeatures);
    free(dataSplits[1].designMatrixY);
    free(dataSplits[1].auxiliaryFeatures);
    free(dataSplits);
    return NULL;
  }
  // lambda in parent node
  double *lambda = leafOutputFunction(
      nrows,
      ncolsDesign,
      ncolsAuxiliary,
      designMatrixY,
      auxiliaryFeatures,
      unitsOfCPIUs,
      treeId,
      nUnits,
      lenOutput);

  // get design matrix and auxiliary features for left and right nodes
  double **designMatrixYLeft = dataSplits[0].designMatrixY;
  double **designMatrixYRight = dataSplits[1].designMatrixY;
  double **auxiliaryFeaturesLeft = dataSplits[0].auxiliaryFeatures;
  double **auxiliaryFeaturesRight = dataSplits[1].auxiliaryFeatures;

  // lambda in left node
  double *lambdaLeft = leafOutputFunction(
      nLeft,
      ncolsDesign,
      ncolsAuxiliary,
      designMatrixYLeft,
      auxiliaryFeaturesLeft,
      unitsOfCPIUs,
      treeId,
      nUnits,
      lenOutput);

  // lambda in right node
  double *lambdaRight = leafOutputFunction(
      nRight,
      ncolsDesign,
      ncolsAuxiliary,
      designMatrixYRight,
      auxiliaryFeaturesRight,
      unitsOfCPIUs,
      treeId,
      nUnits,
      lenOutput);

  // calculate the pseudo likelihood
  double stat = 0.0, statLeft = 0.0, statRight = 0.0;

  for (size_t i = 0; i < nrows; i++)
  {

    stat += designMatrixY[i][ncolsDesign] * log(lambda[0] * auxiliaryFeatures[i][ncolsAuxiliary - 1]) - lambda[0] * auxiliaryFeatures[i][ncolsAuxiliary - 1] - lgamma(designMatrixY[i][ncolsDesign] + 1.0);
    if (designMatrixY[i][splitIndex] <= splitValue)
    {
      statLeft += designMatrixY[i][ncolsDesign] * log(lambdaLeft[0] * auxiliaryFeatures[i][ncolsAuxiliary - 1]) - lambdaLeft[0] * auxiliaryFeatures[i][ncolsAuxiliary - 1] - lgamma(designMatrixY[i][ncolsDesign] + 1.0);
    }
    else
    {
      statRight += designMatrixY[i][ncolsDesign] * log(lambdaRight[0] * auxiliaryFeatures[i][ncolsAuxiliary - 1]) - lambdaRight[0] * auxiliaryFeatures[i][ncolsAuxiliary - 1] - lgamma(designMatrixY[i][ncolsDesign] + 1.0);
    }
  }

  // free memory
  free(lambda);
  free(lambdaLeft);
  free(lambdaRight);
  free(dataSplits[0].designMatrixY);
  free(dataSplits[0].auxiliaryFeatures);
  free(dataSplits[1].designMatrixY);
  free(dataSplits[1].auxiliaryFeatures);
  free(dataSplits);
  // return the pseudo likelihood
  double *out = (double *)calloc(4, sizeof(double));
  out[0] = (statLeft + statRight - stat) / fabs(stat);
  out[1] = statLeft + statRight - stat;
  out[2] = statLeft;
  out[3] = statRight;

  return out;
}

// generic quasi-poisson likelihood function
// Quasi likelihood with interval consideration
double *QuasiPoissonLikelihood(
    fpLeafOutput leafOutputFunction,
    unsigned int splitIndex,
    double splitValue,
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIUs,
    long treeId,
    size_t nUnits,
    size_t lenOutput,
    size_t minNodeSize)
{
  /****************************************************
   *                 Split the Dataset                *
   ****************************************************/
  DecisionTreeData *dataSplits = SplitDataset(
      splitIndex,
      splitValue,
      nrows,
      designMatrixY,
      auxiliaryFeatures);
  // get the number of rows in left and right nodes
  size_t nLeft = dataSplits[0].nrows;
  size_t nRight = dataSplits[1].nrows;

  // check if the number of rows in left and right nodes are zero
  if (nLeft <= minNodeSize || nRight <= minNodeSize)
  {
    free(dataSplits[0].designMatrixY);
    free(dataSplits[0].auxiliaryFeatures);
    free(dataSplits[1].designMatrixY);
    free(dataSplits[1].auxiliaryFeatures);
    free(dataSplits);
    return NULL;
  }
  /******************For dynamic pseudo-risk time estimation**************************/
  // pseudo risk time for each patient each interval in parent, left, right nodes
  double ***pseudoRiskTime = Allocate3DArray(3, nrows, nUnits);

  if (_dynamicRisk == 1)
  {
    if (_verbose >= 2)
    {
      printf("Dynamic pseudo-risk time estimation is enabled\n");
    }
    /******************For dynamic pseudo-risk time estimation**************************/
    // pseudo risk time at current root, left and right nodes
    // pseudo risk time at tree-level is handle in GrowTree() of tree.c
    // no pseudo risk time and population level pseudo risk time are handled by the user input
    double *X = GetCol(auxiliaryFeatures, nrows, ncolsAuxiliary, ncolsAuxiliary - 2);
    double *status = GetCol(auxiliaryFeatures, nrows, ncolsAuxiliary, ncolsAuxiliary - 1);
    double *statusInverse = (double *)malloc(nrows * sizeof(double));

    double *XLeft = GetCol(dataSplits[0].auxiliaryFeatures, nLeft, ncolsAuxiliary, ncolsAuxiliary - 2);
    double *statusLeft = GetCol(dataSplits[0].auxiliaryFeatures, nLeft, ncolsAuxiliary, ncolsAuxiliary - 1);
    double *statusInverseLeft = (double *)malloc(nLeft * sizeof(double));

    double *XRight = GetCol(dataSplits[1].auxiliaryFeatures, nRight, ncolsAuxiliary, ncolsAuxiliary - 2);
    double *statusRight = GetCol(dataSplits[1].auxiliaryFeatures, nRight, ncolsAuxiliary, ncolsAuxiliary - 1);
    double *statusInverseRight = (double *)malloc(nRight * sizeof(double));

    if (statusInverse == NULL || statusInverseLeft == NULL || statusInverseRight == NULL)
    {
      PRINT_LOCATION();
      printf("Memory allocation failed\n");
      exit(1);
    }

    for (int i = 0; i < nrows; i++)
    {
      statusInverse[i] = fabs(status[i]) > 0.99 ? 1.0 : 0.0;
    }

    for (int i = 0; i < nLeft; i++)
    {
      statusInverseLeft[i] = fabs(statusLeft[i]) > 0.99 ? 1.0 : 0.0;
    }

    for (int i = 0; i < nRight; i++)
    {
      statusInverseRight[i] = fabs(statusRight[i]) > 0.99 ? 1.0 : 0.0;
    }

    // calculate the breakpoint of intervals
    double *breakpoints = (double *)malloc((nUnits + 1) * sizeof(double));
    breakpoints[0] = 0;

    for (int i = 0; i < nUnits; i++)
    {
      breakpoints[i + 1] = breakpoints[i] + unitsOfCPIUs[i];
    }

    // fit kaplan-meier curve w.r.t censoring
    KMResult *Gt = KaplanMeier(X, statusInverse, nrows);
    KMResult *GtLeft = KaplanMeier(XLeft, statusInverseLeft, nLeft);
    KMResult *GtRight = KaplanMeier(XRight, statusInverseRight, nRight);
    if (_verbose >= 3)
    {
      printf("Gt:\n");
      PrintArrayDouble(Gt->uniqueTime, Gt->nTime);
      PrintArrayDouble(Gt->survivalProb, Gt->nTime);
      printf("GtLeft:\n");
      PrintArrayDouble(GtLeft->uniqueTime, GtLeft->nTime);
      PrintArrayDouble(GtLeft->survivalProb, GtLeft->nTime);
      printf("GtRight:\n");
      PrintArrayDouble(GtRight->uniqueTime, GtRight->nTime);
      PrintArrayDouble(GtRight->survivalProb, GtRight->nTime);
    }

    free(statusInverse);
    free(statusInverseLeft);
    free(statusInverseRight);

    // placeholder for inverse probability of censoring weight w.r.t time of interest
    KMResult *wtPatient, *wtPatientLeft, *wtPatientRight;
    for (int i = 0; i < nrows; i++)
    {
      wtPatient = Wt(Gt, X[i], status[i]);
      for (int j = 0; j < nUnits; j++)
      {
        pseudoRiskTime[0][i][j] = PseudoRiskTime(wtPatient, breakpoints[j], breakpoints[j + 1]);
      }
      FreeKMResult(wtPatient);
    }
    FreeKMResult(Gt);

    for (int i = 0; i < nLeft; i++)
    {
      wtPatientLeft = Wt(GtLeft, XLeft[i], statusLeft[i]);
      for (int j = 0; j < nUnits; j++)
      {
        pseudoRiskTime[1][i][j] = PseudoRiskTime(wtPatientLeft, breakpoints[j], breakpoints[j + 1]);
      }
      FreeKMResult(wtPatientLeft);
    }
    FreeKMResult(GtLeft);

    for (int i = 0; i < nRight; i++)
    {
      wtPatientRight = Wt(GtRight, XRight[i], statusRight[i]);
      for (int j = 0; j < nUnits; j++)
      {
        pseudoRiskTime[2][i][j] = PseudoRiskTime(wtPatientRight, breakpoints[j], breakpoints[j + 1]);
      }
      FreeKMResult(wtPatientRight);
    }
    FreeKMResult(GtRight);

    if (_verbose >= 3)
    {
      printf("pseudotime parent:\n");
      for (int i = 0; i < nrows; i++)
      {
        printf("i = %d, X = %.3f, Status = %.3f: ", i, X[i], status[i]);
        PrintArrayDouble(pseudoRiskTime[0][i], nUnits);
      }
      printf("pseudotime left:\n");
      for (int i = 0; i < nLeft; i++)
      {
        printf("i = %d, X = %.3f, Status = %.3f: ", i, XLeft[i], statusLeft[i]);
        PrintArrayDouble(pseudoRiskTime[1][i], nUnits);
      }
      printf("pseudotime right:\n");
      for (int i = 0; i < nRight; i++)
      {
        printf("i = %d, X = %.3f, Status = %.3f: ", i, XRight[i], statusRight[i]);
        PrintArrayDouble(pseudoRiskTime[2][i], nUnits);
      }
    }

    free(X);
    free(XLeft);
    free(XRight);
    free(status);
    free(statusLeft);
    free(statusRight);
    free(breakpoints);
  }
  else
  {
    /**************************Inherite the pseudo-risk time from auxiliary features**************************************/
    for (int i = 0; i < nrows; i++)
    {
      for (int j = 0; j < nUnits; j++)
      {
        pseudoRiskTime[0][i][j] = auxiliaryFeatures[i][j + 1];
      }
    }

    for (int i = 0; i < nLeft; i++)
    {
      for (int j = 0; j < nUnits; j++)
      {
        pseudoRiskTime[1][i][j] = dataSplits[0].auxiliaryFeatures[i][j + 1];
      }
    }

    for (int i = 0; i < nRight; i++)
    {
      for (int j = 0; j < nUnits; j++)
      {
        pseudoRiskTime[2][i][j] = dataSplits[1].auxiliaryFeatures[i][j + 1];
      }
    }
  }

  double **lambdas = Allocate2DArray(3, lenOutput);      // lambdas in parent, left, right nodes
  double **Ys = Allocate2DArray(3, lenOutput);           // number of composite events in parent, left, right nodes
  double **rts = Allocate2DArray(3, lenOutput);          // pseudo risk time in parent, left, right nodes
  double **phi = Allocate2DArray(3, nUnits);             // phi in parent, left, right nodes
  double *YTotal = (double *)calloc(3, sizeof(double));  // total number of composite events in parent, left, right nodes
  double *rtTotal = (double *)calloc(3, sizeof(double)); // total pseudo risk time in parent, left, right nodes
  double *betas = (double *)calloc(3, sizeof(double));   // beta in parent, left, right nodes

  /******************For dynamic phi estimation**************************/
  double ***YY = NULL;
  int **lenYY = NULL;

  if (_dynamicPhi == 1)
  {
    YY = Allocate3DArray(3, lenOutput, nrows); // YY = Y/sqrt(rt) for phi estimation in parent, left, right nodes
    lenYY = Allocate2DArrayInt(3, lenOutput);  // lenth of after NA remove for each row of YY for parent, left, right nodes
  }
  /********** calculate interval specific parameters***************/
  for (size_t j = 0; j < lenOutput; j++)
  {
    for (size_t i = 0; i < nrows; i++)
    {
      // make sure the pseudo risk time is not 0 (NA)
      if (fabs(pseudoRiskTime[0][i][j]) > 1e-6)
      {
        Ys[0][j] += designMatrixY[i][ncolsDesign + j];
        rts[0][j] += pseudoRiskTime[0][i][j];
        YTotal[0] += designMatrixY[i][ncolsDesign + j];
        rtTotal[0] += pseudoRiskTime[0][i][j];
        if (_dynamicPhi == 1)
        {
          YY[0][j][lenYY[0][j]++] = designMatrixY[i][ncolsDesign + j] / sqrt(pseudoRiskTime[0][i][j]);
        }
      }
    }
  }

  for (size_t j = 0; j < lenOutput; j++)
  {
    for (size_t i = 0; i < nLeft; i++)
    {
      // make sure the pseudo risk time is not 0 (NA)
      if (fabs(pseudoRiskTime[1][i][j]) > 1e-6)
      {
        Ys[1][j] += dataSplits[0].designMatrixY[i][ncolsDesign + j];
        rts[1][j] += pseudoRiskTime[1][i][j];
        YTotal[1] += dataSplits[0].designMatrixY[i][ncolsDesign + j];
        rtTotal[1] += pseudoRiskTime[1][i][j];
        if (_dynamicPhi == 1)
        {
          YY[1][j][lenYY[1][j]++] = dataSplits[0].designMatrixY[i][ncolsDesign + j] / sqrt(pseudoRiskTime[1][i][j]);
        }
      }
    }
  }

  for (size_t j = 0; j < lenOutput; j++)
  {
    for (size_t i = 0; i < nRight; i++)
    {
      // make sure the pseudo risk time is not 0 (NA)
      if (fabs(pseudoRiskTime[2][i][j]) > 1e-6)
      {
        Ys[2][j] += dataSplits[1].designMatrixY[i][ncolsDesign + j];
        rts[2][j] += pseudoRiskTime[2][i][j];
        YTotal[2] += dataSplits[1].designMatrixY[i][ncolsDesign + j];
        rtTotal[2] += pseudoRiskTime[2][i][j];
        if (_dynamicPhi == 1)
        {
          YY[2][j][lenYY[2][j]++] = dataSplits[1].designMatrixY[i][ncolsDesign + j] / sqrt(pseudoRiskTime[2][i][j]);
        }
      }
    }
  }

  // calculate the beta for parent, left, right nodes
  double alpha = 1.0 / (_k * _k);
  for (size_t i = 0; i < 3; i++)
  {
    betas[i] = alpha / (YTotal[i] / (rtTotal[i] + 1e-9) + 1e-9);
  }

  if (_verbose >= 2)
  {
    printf("beta:\n");
    PrintArrayDouble(betas, 3);
  }
  // calculate the lambda for parent, left, right nodes
  for (size_t i = 0; i < 3; i++)
  {
    if (_verbose >= 2)
    {
      printf("Ys[%ld], rts[%ld]:\n", i, i);
      PrintArrayDouble(Ys[i], lenOutput);
      PrintArrayDouble(rts[i], lenOutput);
    }
    for (size_t j = 0; j < lenOutput; j++)
    {
      lambdas[i][j] = (alpha + Ys[i][j]) / (betas[i] + rts[i][j] + 1e-9);
    }
  }

  if (_verbose >= 2)
  {
    printf("lambda:\n");
    PrintArrayDouble(lambdas[0], lenOutput);
    PrintArrayDouble(lambdas[1], lenOutput);
    PrintArrayDouble(lambdas[2], lenOutput);
  }

  // calculate the phi for parent, left, right nodes, phi is variance of YY * lambda
  for (size_t i = 0; i < 3; i++)
  {
    for (size_t j = 0; j < nUnits; j++)
    {
      if (_dynamicPhi == 1)
      {
        // phi is at the split level
        phi[i][j] = Var(YY[i][j], lenYY[i][j]) / lambdas[i][j];
      }
      else
      {
        // phi is at the tree level(1, forest, tree)
        // note phi at tree and forest level, is implemented at `forest.c` while starting the parallism
        phi[i][j] = _treePhi[treeId][j];
      }
    }
  }

  if (_verbose >= 2)
  {
    printf("phi:\n");
    PrintArrayDouble(phi[0], nUnits);
    PrintArrayDouble(phi[1], nUnits);
    PrintArrayDouble(phi[2], nUnits);
  }
  // calculate the pseudo likelihood
  double stat = 0.0, statLeft = 0.0, statRight = 0.0;

  for (size_t i = 0; i < nrows; i++)
  {
    for (size_t j = 0; j < nUnits; j++)
    {
      // wide format data will include the CPIU with 0 pseudo risk time
      if (fabs(pseudoRiskTime[0][i][j] - 0) > 1e-6)
      {
        stat += 1 / phi[0][j] * (designMatrixY[i][ncolsDesign + j] * log(lambdas[0][j] * pseudoRiskTime[0][i][j]) - lambdas[0][j] * pseudoRiskTime[0][i][j]);
      }
    }
  }

  for (size_t i = 0; i < nLeft; i++)
  {
    for (size_t j = 0; j < nUnits; j++)
    {
      if (fabs(pseudoRiskTime[1][i][j] - 0) > 1e-6)
      {
        statLeft += 1 / phi[1][j] * (dataSplits[0].designMatrixY[i][ncolsDesign + j] * log(lambdas[1][j] * pseudoRiskTime[1][i][j]) - lambdas[1][j] * pseudoRiskTime[1][i][j]);
      }
    }
  }

  for (size_t i = 0; i < nRight; i++)
  {
    for (size_t j = 0; j < nUnits; j++)
    {
      if (fabs(pseudoRiskTime[2][i][j] - 0) > 1e-6)
      {
        statRight += 1 / phi[2][j] * (dataSplits[1].designMatrixY[i][ncolsDesign + j] * log(lambdas[2][j] * pseudoRiskTime[2][i][j]) - lambdas[2][j] * pseudoRiskTime[2][i][j]);
      }
    }
  }

  // free memory
  Free2DArray(lambdas, 3);
  Free2DArray(Ys, 3);
  Free2DArray(rts, 3);
  Free3DArray(YY, 3, lenOutput);
  Free2DArrayInt(lenYY, 3);
  Free2DArray(phi, 3);
  free(YTotal);
  free(rtTotal);
  free(betas);
  free(dataSplits[0].designMatrixY);
  free(dataSplits[0].auxiliaryFeatures);
  free(dataSplits[1].designMatrixY);
  free(dataSplits[1].auxiliaryFeatures);
  free(dataSplits);
  Free3DArray(pseudoRiskTime, 3, nrows);

  // return the quasi likelihood
  double *out = (double *)calloc(4, sizeof(double));
  out[0] = (statLeft + statRight - stat) / fabs(stat);
  out[1] = statLeft + statRight - stat;
  out[2] = statLeft;
  out[3] = statRight;
  if (_verbose >= 2)
  {
    printf("Quasi likelihood: %.3f\n", out[0]);
    PrintArrayDouble(out, 4);
  }
  return out;
}

// Asympotic Test Statistics of Difference of lambdas at each interval
double *AsympoticDiffTestStatsNew(
    fpLeafOutput leafOutputFunction,
    unsigned int splitIndex,
    double splitValue,
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIUs,
    long treeId,
    size_t nUnits,
    size_t lenOutput,
    size_t minNodeSize)
{
  /****************************************************
   *                 Split the Dataset                *
   ****************************************************/
  DecisionTreeData *dataSplits = SplitDataset(
      splitIndex,
      splitValue,
      nrows,
      designMatrixY,
      auxiliaryFeatures);
  // get the number of rows in left and right nodes
  size_t nLeft = dataSplits[0].nrows;
  size_t nRight = dataSplits[1].nrows;

  // check if the number of rows in left and right nodes are zero
  if (nLeft <= minNodeSize || nRight <= minNodeSize)
  {
    free(dataSplits[0].designMatrixY);
    free(dataSplits[0].auxiliaryFeatures);
    free(dataSplits[1].designMatrixY);
    free(dataSplits[1].auxiliaryFeatures);
    free(dataSplits);
    return NULL;
  }

  double **lambdas = Allocate2DArray(3, lenOutput);      // lambdas in parent, left, right nodes
  double **Ys = Allocate2DArray(3, lenOutput);           // number of composite events in parent, left, right nodes
  double **rts = Allocate2DArray(3, lenOutput);          // pseudo risk time in parent, left, right nodes
  double ***YY = Allocate3DArray(3, lenOutput, nrows);   // YY = Y/sqrt(rt) for phi estimation in parent, left, right nodes
  int **lenYY = Allocate2DArrayInt(3, lenOutput);        // lenth of after NA remove for each row of YY for parent, left, right nodes
  double **phi = Allocate2DArray(3, nUnits);             // phi in parent, left, right nodes
  double *YTotal = (double *)calloc(3, sizeof(double));  // total number of composite events in parent, left, right nodes
  double *rtTotal = (double *)calloc(3, sizeof(double)); // total pseudo risk time in parent, left, right nodes
  double *betas = (double *)calloc(3, sizeof(double));   // beta in parent, left, right nodes
  // calculate interval specific lambda

  for (size_t j = 0; j < lenOutput; j++)
  {
    for (size_t i = 0; i < nrows; i++)
    {
      // make sure the pseudo risk time is not 0 (NA)
      if (fabs(auxiliaryFeatures[i][j + 1] - 0) > 1e-6)
      {
        Ys[0][j] += designMatrixY[i][ncolsDesign + j];
        rts[0][j] += auxiliaryFeatures[i][j + 1];
        YTotal[0] += designMatrixY[i][ncolsDesign + j];
        rtTotal[0] += auxiliaryFeatures[i][j + 1];
        YY[0][j][lenYY[0][j]++] = designMatrixY[i][ncolsDesign + j] / sqrt(auxiliaryFeatures[i][j + 1]);
        // check if the ith row is in left node
        if (designMatrixY[i][splitIndex] <= splitValue)
        {
          Ys[1][j] += designMatrixY[i][ncolsDesign + j];
          rts[1][j] += auxiliaryFeatures[i][j + 1];
          YTotal[1] += designMatrixY[i][ncolsDesign + j];
          rtTotal[1] += auxiliaryFeatures[i][j + 1];
          YY[1][j][lenYY[1][j]++] = designMatrixY[i][ncolsDesign + j] / sqrt(auxiliaryFeatures[i][j + 1]);
        }
        else
        {
          Ys[2][j] += designMatrixY[i][ncolsDesign + j];
          rts[2][j] += auxiliaryFeatures[i][j + 1];
          YTotal[2] += designMatrixY[i][ncolsDesign + j];
          rtTotal[2] += auxiliaryFeatures[i][j + 1];
          YY[2][j][lenYY[2][j]++] = designMatrixY[i][ncolsDesign + j] / sqrt(auxiliaryFeatures[i][j + 1]);
        }
      }
    }
  }

  // calculate the beta for parent, left, right nodes
  double alpha = 1.0 / (_k * _k);
  for (size_t i = 0; i < 3; i++)
  {
    betas[i] = alpha / (YTotal[i] / (rtTotal[i] + 1e-9) + 1e-9);
  }

  // calculate the lambda for parent, left, right nodes
  for (size_t i = 0; i < 3; i++)
  {
    for (size_t j = 0; j < lenOutput; j++)
    {
      lambdas[i][j] = (alpha + Ys[i][j]) / (betas[i] + rts[i][j] + 1e-9);
    }
  }

  // calculate the phi for parent, left, right nodes, phi is variance of YY * lambda
  for (size_t i = 0; i < 3; i++)
  {
    for (size_t j = 0; j < nUnits; j++)
    {
      phi[i][j] = Var(YY[i][j], lenYY[i][j]) / lambdas[i][j];
    }
  }

  // calculate the transformed tt for left and right nodes
  double *ttLeft = (double *)calloc(nUnits, sizeof(double));
  double *ttRight = (double *)calloc(nUnits, sizeof(double));
  double *varTTLeft = (double *)calloc(nUnits, sizeof(double));
  double *varTTRight = (double *)calloc(nUnits, sizeof(double));
  for (size_t i = 0; i < nrows; i++)
  {
    for (size_t j = 0; j < nUnits; j++)
    {
      if (fabs(auxiliaryFeatures[i][j + 1] - 0) > 1e-6)
      {
        if (designMatrixY[i][splitIndex] <= splitValue)
        {
          ttLeft[j] += (designMatrixY[i][ncolsDesign + j] - lambdas[0][j] * auxiliaryFeatures[i][j + 1]) / sqrt(auxiliaryFeatures[i][j + 1]);
        }
        else
        {
          ttRight[j] += (designMatrixY[i][ncolsDesign + j] - lambdas[0][j] * auxiliaryFeatures[i][j + 1]) / sqrt(auxiliaryFeatures[i][j + 1]);
        }
      }
    }
  }

  for (size_t j = 0; j < nUnits; j++)
  {
    varTTLeft[j] = phi[1][j] * lambdas[0][j] / sqrt(lenYY[1][j]);
    varTTRight[j] = phi[2][j] * lambdas[0][j] / sqrt(lenYY[2][j]);
    ttLeft[j] = ttLeft[j] / lenYY[1][j];
    ttRight[j] = ttRight[j] / lenYY[2][j];
  }

  // calculate the test statistics
  double ZZ = (Sum(ttLeft, nUnits) - Sum(ttRight, nUnits)) / sqrt(Sum(varTTLeft, nUnits) + Sum(varTTRight, nUnits));
  double temp = 2.0f * gsl_cdf_ugaussian_Q(fabs(ZZ));
  double *out = (double *)calloc(4, sizeof(double));
  out[0] = -1.0 * log10(temp);
  out[1] = fabs(ZZ);
  out[2] = ZZ;
  out[3] = ZZ;
  // free memory
  Free2DArray(lambdas, 3);
  Free2DArray(Ys, 3);
  Free2DArray(rts, 3);
  Free3DArray(YY, 3, lenOutput);
  Free2DArrayInt(lenYY, 3);
  Free2DArray(phi, 3);
  free(YTotal);
  free(rtTotal);
  free(betas);
  free(ttLeft);
  free(ttRight);
  free(varTTLeft);
  free(varTTRight);
  free(dataSplits[0].designMatrixY);
  free(dataSplits[0].auxiliaryFeatures);
  free(dataSplits[1].designMatrixY);
  free(dataSplits[1].auxiliaryFeatures);
  free(dataSplits);

  return out;
}

// Generalized Estimating Equation
double *GEERule(
    fpLeafOutput leafOutputFunction,
    unsigned int splitIndex,
    double splitValue,
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIUs,
    long treeId,
    size_t nUnits,
    size_t lenOutput,
    size_t minNodeSize)
{
  /****************************************************
   *                 Split the Dataset                *
   ****************************************************/
  DecisionTreeData *dataSplits = SplitDataset(
      splitIndex,
      splitValue,
      nrows,
      designMatrixY,
      auxiliaryFeatures);
  // get the number of rows in left and right nodes
  size_t nLeft = dataSplits[0].nrows;
  size_t nRight = dataSplits[1].nrows;

  // check if the number of rows in left and right nodes are zero
  if (nLeft <= minNodeSize || nRight <= minNodeSize)
  {
    free(dataSplits[0].designMatrixY);
    free(dataSplits[0].auxiliaryFeatures);
    free(dataSplits[1].designMatrixY);
    free(dataSplits[1].auxiliaryFeatures);
    free(dataSplits);
    return NULL;
  }
  /******************For dynamic pseudo-risk time estimation**************************/
  // pseudo risk time for each patient each interval in parent, left, right nodes
  double **pseudoRiskTime = Allocate2DArray(nrows, nUnits);
  if (_dynamicRisk == 1)
  {
    if (_verbose >= 2)
    {
      printf("Dynamic pseudo-risk time estimation is enabled\n");
    }
    /******************For dynamic pseudo-risk time estimation**************************/
    // pseudo risk time at current root, left and right nodes
    // pseudo risk time at tree-level is handle in GrowTree() of tree.c
    // no pseudo risk time and population level pseudo risk time are handled by the user input
    double *X = GetCol(auxiliaryFeatures, nrows, ncolsAuxiliary, ncolsAuxiliary - 2);
    double *status = GetCol(auxiliaryFeatures, nrows, ncolsAuxiliary, ncolsAuxiliary - 1);
    double *statusInverse = (double *)malloc(nrows * sizeof(double));

    if (statusInverse == NULL)
    {
      PRINT_LOCATION();
      printf("Memory allocation failed\n");
      exit(1);
    }

    for (int i = 0; i < nrows; i++)
    {
      if (fabs(status[i]) > 0.99)
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
    breakpoints[0] = 0;

    for (int i = 0; i < nUnits; i++)
    {
      breakpoints[i + 1] = breakpoints[i] + unitsOfCPIUs[i];
    }

    // fit kaplan-meier curve w.r.t censoring
    KMResult *Gt = KaplanMeier(X, statusInverse, nrows);
    if (_verbose >= 2)
    {
      printf("Gt:\n");
      PrintArrayDouble(Gt->uniqueTime, Gt->nTime);
      PrintArrayDouble(Gt->survivalProb, Gt->nTime);
    }
    free(statusInverse);

    // placeholder for inverse probability of censoring weight w.r.t time of interest
    KMResult *wtPatient;
    for (int i = 0; i < nrows; i++)
    {
      wtPatient = Wt(Gt, X[i], status[i]);
      for (int j = 0; j < nUnits; j++)
      {
        pseudoRiskTime[i][j] = PseudoRiskTime(wtPatient, breakpoints[j], breakpoints[j + 1]);
      }
      FreeKMResult(wtPatient);
    }
    FreeKMResult(Gt);

    if (_verbose >= 3)
    {
      printf("pseudotime parent:\n");
      for (int i = 0; i < nrows; i++)
      {
        printf("i = %d, X = %.3f, Status = %.3f: ", i, X[i], status[i]);
        PrintArrayDouble(pseudoRiskTime[i], nUnits);
      }
    }

    free(X);
    free(status);
    free(breakpoints);
  }
  else
  {
    /**************************Inherite the pseudo-risk time from auxiliary features**************************************/
    for (int i = 0; i < nrows; i++)
    {
      for (int j = 0; j < nUnits; j++)
      {
        pseudoRiskTime[i][j] = auxiliaryFeatures[i][j + 1];
      }
    }
  }
  /********************GEE modeling*************************/
  // data preparation

  MATRIX *xin, *yin, *idin, *nin, *offsetin, *betain;
  MATRIX *naivvar = NULL, *robvar = NULL, *tempmat1 = NULL, *tempmat2 = NULL;

  xin = _interaction == 1 ? VC_GEE_create_matrix(nrows * nUnits, 2 * nUnits, PERMANENT) : VC_GEE_create_matrix(nrows * nUnits, nUnits + 1, PERMANENT);
  yin = VC_GEE_create_matrix(nrows * nUnits, 1, PERMANENT);
  idin = VC_GEE_create_matrix(nrows * nUnits, 1, PERMANENT);
  nin = VC_GEE_create_matrix(nrows * nUnits, 1, PERMANENT);
  offsetin = VC_GEE_create_matrix(nrows * nUnits, 1, PERMANENT);
  betain = _interaction == 1 ? VC_GEE_create_matrix(2 * nUnits, 1, PERMANENT) : VC_GEE_create_matrix(nUnits + 1, 1, PERMANENT);

  for (int i = 0; i < nrows; i++)
  {
    for (int j = 0; j < (nUnits - 1); j++)
    {
      // intercept [0]
      xin->data[i * nUnits + j][0] = 1;
      // interval [1]: [nUnits - 1]
      xin->data[i * nUnits + j][j + 1] = 1;
      // var [nUnits]
      xin->data[i * nUnits + j][nUnits] = designMatrixY[i][splitIndex] <= splitValue ? 1 : 0;
      // var * interval [Units + 1] : [2 * nUnits -1]
      if (_interaction == 1)
        xin->data[i * nUnits + j][nUnits + j + 1] = xin->data[i * nUnits + j][nUnits];
      yin->data[i * nUnits + j][0] = designMatrixY[i][ncolsDesign + j];
      idin->data[i * nUnits + j][0] = i;
      nin->data[i * nUnits + j][0] = 1;
      offsetin->data[i * nUnits + j][0] = log(pseudoRiskTime[i][j] + 1e-5);
    }
    xin->data[(i + 1) * nUnits - 1][0] = 1;
    xin->data[(i + 1) * nUnits - 1][nUnits] = designMatrixY[i][splitIndex] <= splitValue ? 1 : 0;
    yin->data[(i + 1) * nUnits - 1][0] = designMatrixY[i][ncolsDesign + nUnits - 1];
    idin->data[(i + 1) * nUnits - 1][0] = i;
    nin->data[(i + 1) * nUnits - 1][0] = 1;
    offsetin->data[(i + 1) * nUnits - 1][0] = log(pseudoRiskTime[i][nUnits - 1] + 1e-5);
  }

  int nobs = nrows * nUnits;
  int p = _interaction == 1 ? 2 * nUnits : nUnits + 1;
  int parmvec[3];
  parmvec[0] = 2;
  parmvec[1] = 2;
  parmvec[2] = 1;

  int M_parm = 2;
  int maxsz = nUnits;
  int S_iter = 50;
  int silent = 1;
  int scale_fix = 0;
  int compatflag = 0;

  double tol = 1e-4;
  double S_phi = 1.0;

  MATRIX *Rin = NULL;

  if (_Rin != NULL)
  {
    Rin = VC_GEE_matcopy(&_Rin);
  }

  int errorState = Cgee(
      &xin,
      &yin,
      &idin,
      &nin,
      &offsetin,
      &nobs,
      &p,
      parmvec,
      &M_parm,
      &betain,
      &naivvar,
      &robvar,
      &Rin,
      &S_phi,
      &tol,
      &maxsz,
      &S_iter,
      &silent,
      &scale_fix,
      &compatflag);

  VC_GEE_destroy_matrix(&Rin);

  // return the quasi likelihood

  if (errorState != 0)
  {
    if (_verbose >= 3)
      printf("Warning: GEE fitting failed; Error code: %d\n", errorState);
    VC_GEE_destroy_matrix(&naivvar);
    VC_GEE_destroy_matrix(&robvar);
    VC_GEE_destroy_matrix(&betain);
    VC_GEE_destroy_matrix(&xin);
    VC_GEE_destroy_matrix(&yin);
    VC_GEE_destroy_matrix(&idin);
    VC_GEE_destroy_matrix(&nin);
    VC_GEE_destroy_matrix(&offsetin);
    free(dataSplits[0].designMatrixY);
    free(dataSplits[0].auxiliaryFeatures);
    free(dataSplits[1].designMatrixY);
    free(dataSplits[1].auxiliaryFeatures);
    free(dataSplits);
    Free2DArray(pseudoRiskTime, nrows);
    return NULL;
  }

  double *out = (double *)calloc(4, sizeof(double));

  if (_interaction == 1)
  {
    tempmat1 = VC_GEE_extract_rows(&betain, nUnits, 2 * nUnits - 1);
    tempmat2 = VC_GEE_ident(nUnits);
    MATRIX *A = VC_GEE_matmult(&tempmat2, &tempmat1);
    // MATRIX *A = VC_GEE_matmult(VC_GEE_ident(nUnits), VC_GEE_extract_rows(betain, nUnits, 2 * nUnits - 1));                        // (nUnits * nUnits) %*% (nUnits * 1) = (nUnits * 1)
    tempmat1 = VC_GEE_extract_rows(&robvar, nUnits, 2 * nUnits - 1);
    tempmat2 = VC_GEE_extract_cols(&tempmat1, nUnits, 2 * nUnits - 1);
    MATRIX *B = matrix_inverse(&tempmat2); // (nUnits * nUnits)
    // MATRIX *B = matrix_inverse(VC_GEE_extract_cols(VC_GEE_extract_rows(robvar, nUnits, 2 * nUnits - 1), nUnits, 2 * nUnits - 1)); // (nUnits * nUnits)
    make_permanent(A);
    make_permanent(B);
    tempmat1 = VC_GEE_transp(&A);
    tempmat2 = VC_GEE_matmult(&B, &A);
    MATRIX *OUT = VC_GEE_matmult(&tempmat1, &tempmat2); // (1 * nUnits) %*% (nUnits * nUnits) %*% (nUnits * 1) = (1 * 1)
    // MATRIX *OUT = VC_GEE_matmult(VC_GEE_transp(A), VC_GEE_matmult(B, A)); // (1 * nUnits) %*% (nUnits * nUnits) %*% (nUnits * 1) = (1 * 1)

    double temp = gsl_cdf_chisq_Q(OUT->data[0][0], nUnits);

    out[0] = -1.0 * log10(temp);
    out[1] = -1.0 * log10(temp);
    out[2] = OUT->data[0][0];
    out[3] = temp;

    if (_verbose >= 2)
    {
      printf("Gee output: %.3f\n", out[0]);
      PrintArrayDouble(out, 4);
      printf("nLeft: %ld, nRight: %ld\n", nLeft, nRight);
    }

    if (_verbose >= 3)
    {
      printf("Matrix beta:\n");
      PrintMatrix(&betain);
      printf("Matrix naivvar:\n");
      PrintMatrix(&naivvar);
      printf("Matrix A:\n");
      PrintMatrix(&A);
      printf("Matrix B:\n");
      PrintMatrix(&B);
      printf("Matrix OUT:\n");
      PrintMatrix(&OUT);
    }

    // free memory
    VC_GEE_destroy_matrix(&A);
    VC_GEE_destroy_matrix(&B);
    VC_GEE_destroy_matrix(&OUT);
  }
  else
  {
    double coefBeta = betain->data[nUnits][0];
    double varBeta = robvar->data[nUnits][nUnits];
    double testStat = coefBeta / sqrt(varBeta);
    double pValue = 2.0f * gsl_cdf_ugaussian_Q(fabs(testStat));
    out[0] = -1.0 * log10(pValue);
    out[1] = -1.0 * log10(pValue);
    out[2] = testStat;
    out[3] = pValue;
  }

  VC_GEE_destroy_matrix(&naivvar);
  VC_GEE_destroy_matrix(&robvar);
  VC_GEE_destroy_matrix(&betain);
  VC_GEE_destroy_matrix(&xin);
  VC_GEE_destroy_matrix(&yin);
  VC_GEE_destroy_matrix(&idin);
  VC_GEE_destroy_matrix(&nin);
  VC_GEE_destroy_matrix(&offsetin);
  free(dataSplits[0].designMatrixY);
  free(dataSplits[0].auxiliaryFeatures);
  free(dataSplits[1].designMatrixY);
  free(dataSplits[1].auxiliaryFeatures);
  free(dataSplits);
  Free2DArray(pseudoRiskTime, nrows);
  return out;
}

// get the best split
SplitPoints *FindBestSplit(
    fpSplitFunction splitFunction,
    fpLeafOutput leafOutputFunction,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIU,
    long treeId,
    size_t nUnits,
    size_t lenOutput,
    size_t mtry,
    size_t nsplits,
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    size_t minNodeSize,
    unsigned int seed)
{
  int bestSplitIndex = -1;
  double bestSplitValue = 0.0;
  double *bestSplitStat = calloc(4, sizeof(double));

  // randomly select non-repeated mtry variables
  int *varIndexs = GetSeqInt(0, (int)ncolsDesign - 1, 1);
  if (_verbose >= 2)
  {
    printf("ncolDesign: %ld\n", ncolsDesign);
    printf("varIndexs:\n");
  }
  int *varIndexSample = SampleInt(varIndexs, ncolsDesign, mtry, 0, seed);
  free(varIndexs);

  if (_verbose >= 2)
  {
    printf("varIndexSample:\n");
    PrintArrayInt(varIndexSample, mtry);
  }

  // randomly select non-repeated nsplits values for each variable
  ElementStruct **splitsOfmtry = (ElementStruct **)malloc(mtry * sizeof(ElementStruct *));

  size_t maxSplits = 0;
  for (size_t i = 0; i < mtry; i++)
  {
    int varIndex = varIndexSample[i];
    double *varValues = GetCol(designMatrixY, nrows, ncolsDesign, varIndex);
    ElementStruct *varValuesNaRemoved = RemoveNA(varValues, nrows);
    splitsOfmtry[i] = SplitCandidates(
        varValuesNaRemoved->Elements,
        varValuesNaRemoved->nElements,
        nsplits);
    maxSplits += splitsOfmtry[i]->nElements;
    free(varValues);
    free(varValuesNaRemoved->Elements);
    free(varValuesNaRemoved);
  }

  double *pValues = (double *)calloc(maxSplits, sizeof(double));

  if (_verbose >= 2)
  {
    printf("mtry: %ld\n", mtry);
    for (size_t i = 0; i < mtry; i++)
    {
      printf("varIndexSample[%ld]: %d\n", i, varIndexSample[i]);
      printf("splitsOfmtry[%ld]:\n", i);
      PrintArrayDouble(splitsOfmtry[i]->Elements, splitsOfmtry[i]->nElements);
    }
  }

  int totalValidSplits = 0;

  // find the best split
  for (size_t i = 0; i < mtry; i++)
  {
    for (size_t j = 0; j < splitsOfmtry[i]->nElements; j++)
    {
      // calculate the split stat
      double *splitStat = splitFunction(
          leafOutputFunction,
          varIndexSample[i],
          splitsOfmtry[i]->Elements[j],
          nrows,
          ncolsDesign,
          ncolsAuxiliary,
          designMatrixY,
          auxiliaryFeatures,
          unitsOfCPIU,
          treeId,
          nUnits,
          lenOutput,
          minNodeSize);

      if (splitStat != NULL)
      {
        pValues[totalValidSplits++] = splitStat[3];
        if (splitStat[0] > bestSplitStat[0])
        {
          bestSplitIndex = varIndexSample[i];
          bestSplitValue = splitsOfmtry[i]->Elements[j];
          memcpy(bestSplitStat, splitStat, 4 * sizeof(double));
        }
      }
      free(splitStat);
    }
  }

  SplitPoints *bestSplit = (SplitPoints *)malloc(sizeof(SplitPoints));
  bestSplit->splitIndex = bestSplitIndex;
  bestSplit->splitValue = bestSplitValue;
  bestSplit->splitStat = bestSplitStat;

  pValues = realloc(pValues, totalValidSplits * sizeof(double));

  if (_gee == 1 && totalValidSplits > 1 && strcmp(_padjust, "none") != 0)
  {
    double *pValueAdj = PAdjust(pValues, totalValidSplits, _padjust);
    if (pValueAdj != NULL)
    {
      double minPvalueAdj = Min(pValueAdj, totalValidSplits);
      double negLog10minPvalueAdj = -1.0 * log10(minPvalueAdj);
      bestSplit->splitStat[0] = negLog10minPvalueAdj;
      free(pValueAdj);
    }
  }

  // free the memory
  free(varIndexSample);
  for (size_t i = 0; i < mtry; i++)
  {
    free(splitsOfmtry[i]->Elements);
    free(splitsOfmtry[i]);
  }
  free(splitsOfmtry);
  free(pValues);

  return bestSplit;
}
