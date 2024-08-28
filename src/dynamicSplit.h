/*
@author Yu Wang
*/

#ifndef SPLIT_H
#define SPLIT_H

#include "utils.h"
#include "surv.h"

// typedef function pointer for code readability
typedef double *(*fpLeafOutput)(
    size_t,    // nrow
    size_t,    // ncolsDesign
    size_t,    // ncolsAuxiliary
    double **, // designMatrixY
    double **, // auxiliaryFeatures
    double *,  // unitsOfCPIU
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
    double *,     // unitsOfCPIUS
    size_t,       // nUnits
    size_t        // lenOutput
);


// struct to store the data for one node
typedef struct
{
    size_t nrows;
    double **designMatrixY;
    double **auxiliaryFeatures;
} DecisionTreeData;

// split a dataset into two parts according to splitIndex and splitValue
DecisionTreeData *SplitDataset(
    unsigned varIndex,                 // var index
    double splitValue,              // best cutoff
    size_t nrows,                   // nrows of design matrix, number of variables
    double **designMatrixY,   // design matrix (X) and response Y
    double **auxiliaryFeatures // auxiliaryFeatures
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
    size_t nUnits,
    size_t lenOutput,
    size_t mtry,
    size_t nsplits,
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    unsigned int seed);



// output of a leaf node
double *LeafOutput(
    size_t nrow,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIU,
    size_t nUnits,
    size_t lenOutput);

double *LeafOutputInterval(
    size_t nrow,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIU,
    size_t nUnits,
    size_t lenOutput);

    
// compute the pseudo-poisson likelihood based on pool estimation
double *PseudoPoissonLikelihood(
    fpLeafOutput leafOutputFunction,
    unsigned int splitIndex,
    double splitValue,
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIU,
    size_t nUnits,
    size_t lenOutput
);

// compute the pseudo-poisson likelihood based on interval estimation 
double *PseudoPoissonLikelihoodInterval(
    fpLeafOutput leafOutputFunction,
    unsigned int splitIndex,
    double splitValue,
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double **designMatrixY,
    double **auxiliaryFeatures,
    double *unitsOfCPIUS,
    size_t nUnits,
    size_t lenOutput
);

#endif