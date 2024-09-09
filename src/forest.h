/*
@author Yu Wang
*/

#ifndef FOREST_H
#define FOREST_H

#include "utils.h"
#include <omp.h>
#include "tree.h"
#include <dirent.h>   // for mkdir
#include <sys/stat.h> // for mkdir
#include <unistd.h>   // for access

extern int _nPerms;

typedef struct RandomSurvivalForest
{
    size_t nrowsDesign;
    size_t ncolsDesign;
    unsigned int *varIDs;
    size_t nVars;
    size_t nTrees;
    double *unitsOfCPIU;
    size_t nUnits;
    size_t lenOutput;
    size_t maxDepth;
    size_t minNodeSize;
    double minGain;
    size_t mtry;
    size_t nsplits;
    unsigned int seed;

    double **predicted;
    double **oobPredicted;

    double oobMSE;
    double inbagMSE;

    double *vimpStat;      // same lenght as ncolDesign
    double *vimpFreq;      // same length as ncolDesign
    double *likelihoodsum; // quasi-likelihood sum for each tree
    double **vimpPermuted; //permutation importance, the difference between the original and permuted quasi-likelihood sum, same colnumns numbers as nVars;

    int **bagMatrix;

    DecisionTreeNode **forest;

} RandomSurvivalForest;

// forest
RandomSurvivalForest *RandomForest(
    fpSplitFunction splitFunction,
    fpLeafOutput leafOutputFunction,
    double **designMatrixY,
    double **auxiliaryFeatures,
    size_t nrow,
    size_t ncolsDesign,    // there are existing dummy variable coding for categorial variables, use varIDs to recover the vimpPermuted
    size_t ncolsAuxiliary,
    unsigned int  *varIDs, // a vector of variable IDs, from 0 to (nVars -1). the consecutive repeated Ids indicates a single categorial variable when calculating vimpPermuted
    size_t nVars,
    double *unitsOfCPIU,
    size_t nUnits,    // length of unitsOfCPIU
    size_t lenOutput, // length of output at leaf node
    size_t maxDepth,
    size_t minNodeSize,
    double minGain,
    size_t mtry,
    size_t nsplits,
    size_t nTrees,
    unsigned int seed);

double **InternalForestOOBPredict(
    DecisionTreeNode **forest,
    size_t nTrees,
    double **designMatrix,
    size_t nrowsDesign,
    size_t lenOutput,
    int **bagMatrix);

// prediction from forest
double **ForestPredict(
    DecisionTreeNode **forest,
    size_t nTrees,
    double **designMatrix,
    size_t nrowsDesign,
    size_t lenOutput);

double **SurvivalForestPredict(
    RandomSurvivalForest *survivalForest,
    double **designMatrix,
    size_t nrowsDesign);

// gather vimp vector for a forest
void VimpForest(
    DecisionTreeNode **forest,
    size_t nTrees,
    double **vimpStat,
    double **vimpFreq);
// free forest
void FreeSurvivalForest(RandomSurvivalForest *forest);

// save forest
void SaveForest(
    DecisionTreeNode **forest,
    size_t nTrees,
    size_t lenOutput,
    char *path);

// save survival forest
void SaveSurvivalForest(
    RandomSurvivalForest *forest,
    char *path);

// load forest
DecisionTreeNode **LoadForest(
    char *path,
    size_t nTrees);

// load survival forest
RandomSurvivalForest *LoadSurvivalForest(char *path);

// print all tree paths of a forest
void PrintForestPaths(RandomSurvivalForest *forest, FILE *file); 

#endif
