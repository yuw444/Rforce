/*
@author Yu Wang
*/

#ifndef SPLIT_H
#define SPLIT_H

#include "utils.h"
#include "surv.h"
#include "gee.h"
#include "pvalue.h"

extern int _k;

//default: calculate for each tree
extern int _noPseudo;    // make sure the auxiliary input no pseudo-risk
extern int _pseudoRisk1; // use the pseudo-risk time estimate for the entire population
extern int _pseudoRisk2; // use the pseudo-risk time estimate at each tree
extern int _dynamicRisk; // pseudo-risk time estimate at each split

//default: calculate for each tree
extern int _noPhi;       // phi = 1 all the time
extern int _phi1;        // phi is estimated for the entire population
extern int _phi2;        // phi is estimated at each tree
extern int _dynamicPhi;  // phi is calculated at each split

extern int _longformat;   // the original RF-SLAM
extern int _gee;          // the proposed GEE
extern int _interaction;  // interaction term
extern int _asympotic;    // TBD asympotic test

extern int _verbose;      // verbose for debugging; >=1 is on; 0 is off

extern double **_treePhi;  // store the phi for each tree
extern MATRIX *_Rin;       // store the correlation matrix among the CPIU

extern char *_padjust;    // p-value adjustment method

// typedef function pointer for code readability
typedef double *(*fpLeafOutput)(
    size_t,    // nrow
    size_t,    // ncolsDesign
    size_t,    // ncolsAuxiliary
    double **, // designMatrixY
    double **, // auxiliaryFeatures
    double *,  // unitsOfCPIU
    long,      // treeId
    size_t,    // nUnits
    size_t     // lenOutput
);

typedef double *(*fpSplitFunction)(
    fpLeafOutput, // leafOutputFunction
    unsigned int, // splitIndex
    double,       // splitValue
    size_t,       // nrows
    size_t,       // ncolsDesign
    size_t,       // ncolsAuxiliary
    double **,    // designMatrixY
    double **,    // auxiliaryFeatures
    double *,     // unitsOfCPIUs
    long,         // treeId
    size_t,       // nUnits
    size_t,       // lenOutput
    size_t        // minNodeSize
);

// struct to store the data for one node
typedef struct
{
    size_t nrows;
    double **designMatrixY;
    double **auxiliaryFeatures;
} DecisionTreeData;

// find the left and right node size after split
size_t *SplitSizes(
    unsigned int varIndex,     // var index
    double splitValue,         // best cutoff
    size_t nrows,              // nrows of design matrix, number of variables
    double **designMatrixY,    // design matrix (X) and response Y
    double **auxiliaryFeatures // auxiliaryFeatures
);

// split a dataset into two parts according to splitIndex and splitValue
DecisionTreeData *SplitDataset(
    unsigned varIndex,         // var index
    double splitValue,         // best cutoff
    size_t nrows,              // nrows of design matrix, number of variables
    double **designMatrixY,    // design matrix (X) and response Y
    double **auxiliaryFeatures // auxiliaryFeatures
);

// get bootStrapMatrix
int **BagMatrix(
    double *patientIds,
    size_t nrowsDesign,
    size_t nTrees,
    unsigned int seed);

// get bootStrapSample
DecisionTreeData *BootStrapSample(
    double **designMatrixY,
    double **auxiliaryFeatures,
    size_t nrowsDesign,
    int *bootVector // length nrowsDesign
);

// struct to store the best split
typedef struct
{
    int splitIndex;
    double splitValue;
    double *splitStat;
} SplitPoints;

// find the best split for a node
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
    unsigned int seed);

// using the root node pseudotime to estimate pooled lambdas, no interval, for original RF-SLAM only
// with consideration of bayesian estimator, alpha, and beta
double *LeafOutput(
    size_t nrow,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIU,
    long treeId,
    size_t nUnits,
    size_t lenOutput);

// using the root node pseudotime to estimate lambdas at each interval
// with consideration of bayesian estimator, alpha, and beta
double *LeafOutputInterval(
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIUs,
    long treeId,
    size_t nUnits,
    size_t lenOutput);

// compute the Poisson likelihood based on pool estimation, for original RF-SLAM only
// return: percentage gain in the likelihood; parent, left and right likelihood
double *PoissonLikelihood(
    fpLeafOutput leafOutputFunction,
    unsigned int splitIndex,
    double splitValue,
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIU,
    long treeId,
    size_t nUnits,
    size_t lenOutput,
    size_t minNodeSize);

// generic quasi-poisson likelihood function
// return the percentage gain in the likelihood; gain, left and right likelihood
double *QuasiPoissonLikelihood(
    fpLeafOutput leafOutputFunction,
    unsigned int splitIndex,
    double splitValue,
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIU,
    long treeId,
    size_t nUnits,
    size_t lenOutput,
    size_t minNodeSize);

// Asympotic Test Statistics of Difference of lambdas at each interval
// return the test statistics -log10(p-value); fabs(ZZ); p-value; p-value 
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
    size_t minNodeSize);

// Wald Test Statistics for testing GEE coefficients
// return -log10(p-value); Wald test statistics; p-value; p-value
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
    size_t minNodeSize);

enum PhiType
{
    phiNone,   // phi = 1
    phiForest, // calculate for the entire forest
    phiTree,   // calculate at each tree
    phiSplit   // calcualte at each split
};

enum PseudoRiskType
{
    riskNone,   // no risk
    riskForest, // calculate for the entire forest
    riskTree,   // calculate at each tree
    riskSplit   // calculate at each split
};

#endif