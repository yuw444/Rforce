#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include "convert.h"
#include "forest.h"
#include <omp.h>

enum SplitFunctionIndex
{
    QuasiPoissonLikelihoodIndex = 0,
    GEERuleIndex = 1,
    PoissonLikelihoodIndex = 2
};

size_t CountNodes(DecisionTreeNode *root)
{
    if (root == NULL)
    {
        return 0;
    }
    else
    {
        return 1 + CountNodes(root->leftChild) + CountNodes(root->rightChild);
    }
}

void SaveTreeArray(
    DecisionTreeNode *root,
    size_t *rowIndex,
    double ***forestMatrix // nodesCounts * colCounts * nTrees
)
{
    if (root == NULL)
    {
        return;
    }
    else
    {
        size_t treeId = root->treeId;
        size_t nodeId = root->nodeId;
        int leftDaughterNode = root->leftChild == NULL ? -1 : root->leftChild->nodeId;
        int rightDaughterNode = root->rightChild == NULL ? -1 : root->rightChild->nodeId;

        // save the treeId, nodeId, flag, leftDaughter, rightDaughter, splitIndex, splitValue, splitStat, sizeLR, output
        (*forestMatrix)[*rowIndex][0] = treeId * 1.0f;
        (*forestMatrix)[*rowIndex][1] = nodeId * 1.0f;
        (*forestMatrix)[*rowIndex][2] = root->flag * 1.0f;
        (*forestMatrix)[*rowIndex][3] = root->sizeLR[2] * 1.0f;
        (*forestMatrix)[*rowIndex][4] = leftDaughterNode * 1.0f;
        (*forestMatrix)[*rowIndex][5] = rightDaughterNode * 1.0f;
        (*forestMatrix)[*rowIndex][6] = root->splitIndex * 1.0f;
        (*forestMatrix)[*rowIndex][7] = root->splitValue * 1.0f;
        for (size_t i = 0; i < root->lenOutput; i++)
        {
            (*forestMatrix)[*rowIndex][8 + i] = root->output[i];
        }
        (*rowIndex)++;
        SaveTreeArray(root->leftChild, rowIndex, forestMatrix);
        SaveTreeArray(root->rightChild, rowIndex, forestMatrix);
    }
}

SEXP R_Cummax(SEXP x)
{
    double *xptr = REAL(x);
    size_t n = Rf_length(x);
    double *result = Cummax(xptr, n);
    SEXP res = DoublePtrToRVector(result, n);
    free(result);
    return res;
}

SEXP R_ColsPermute(SEXP x, SEXP colsToPermute, SEXP seed)
{

    int nrows = Rf_nrows(x);
    int ncols = Rf_ncols(x);
    // Rprintf("nrows: %d, ncols: %d\n", nrows, ncols);
    int ncolsToPermute = Rf_length(colsToPermute);
    double **xptr = RMatrixToDoublePtr(x);
    // for(int i = 0; i < nrows; i++)
    // {
    //     for(int j = 0; j < ncols; j++)
    //     {
    //         Rprintf("%f ", xptr[i][j]);
    //     }
    //     Rprintf("\n");
    // }
    unsigned int *colsToPermutePtr = (unsigned int *)INTEGER(colsToPermute);
    int seedInt = INTEGER(seed)[0];

    double **result = ColsPermute(xptr, nrows, ncols, colsToPermutePtr, ncolsToPermute, seedInt);
    // for(int i = 0; i < nrows; i++)
    // {
    //     for(int j = 0; j < ncols; j++)
    //     {
    //         Rprintf("%f ", result[i][j]);
    //     }
    //     Rprintf("\n");
    // }
    SEXP res = DoublePtrToRMatrix(result, nrows, ncols);
    return res;
}

// try openmp backend with a simple parallel for loop in C with reduction(+:sum)
// #pragma omp parallel for reduction(+:sum)
SEXP R_Sum(SEXP x, SEXP nthreads)
{
    double *xptr = REAL(x);
    size_t n = Rf_length(x);
    int nthreads0 = INTEGER(nthreads)[0];
    double sum = 0;
    omp_set_num_threads(nthreads0);
    // int nthreads1 = omp_get_max_threads();
#pragma omp parallel for reduction(+ : sum)
    for (size_t i = 0; i < n; i++)
    {
        sum += xptr[i];
        printf("thread id: %d/%d, i: %ld, sum: %f\n", omp_get_thread_num(), omp_get_num_threads(), i, sum);
    }
    // Rprintf("number of thread allocated: %d, requested: %d\n", nthreads1, nthreads0);
    return Rf_ScalarReal(sum);
}

// matrix add
SEXP R_MatrixAdd(SEXP x, SEXP y)
{
    size_t n = Rf_length(y);
    int nrow = Rf_nrows(x);
    int ncol = Rf_ncols(x);
    printf("nrow: %d, ncol: %d\n", nrow, ncol);
    SEXP rst = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
    double *xptr = REAL(x);
    double *yptr = REAL(y);
    double *rstptr = REAL(rst);
    omp_set_num_threads(4);

    // #pragma omp parallel
    // {
    //     printf("thread %d/%d\n", omp_get_thread_num(), omp_get_max_threads());
    // }

#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        rstptr[i] = xptr[i] + yptr[i];
        printf("thread id: %d, i: %ld, sum: %f\n", omp_get_thread_num(), i, rstptr[i]);
    }
    UNPROTECT(1);

    // #pragma omp parallel num_threads(2) // Set 2 threads for this parallel region
    // {
    //     int thread_id = omp_get_thread_num();
    //     // printf("Thread ID in parallel region with 2 threads: %d\n", thread_id);
    // }

    Rprintf("number of thread requested/allocated: %d/%d\n", 4L, omp_get_max_threads());
    return rst;
}

// return a list of vectors to R

SEXP R_ListOfVectors()
{
    SEXP list = PROTECT(Rf_allocVector(VECSXP, 3));

    double v1ptr[3];
    int v2ptr[4];
    char **v3ptr;
    v3ptr = (char **)malloc(5 * sizeof(char *));
    for (size_t i = 0; i < 5; i++)
    {
        v3ptr[i] = (char *)malloc(10 * sizeof(char));
    }

    for (size_t i = 0; i < 3; i++)
    {
        v1ptr[i] = i * 1.0f;
    }

    for (size_t i = 0; i < 4; i++)
    {
        v2ptr[i] = i;
    }

    for (size_t i = 0; i < 5; i++)
    {
        sprintf(v3ptr[i], "string %ld", i);
    }

    SET_VECTOR_ELT(list, 0, DoublePtrToRVector(v1ptr, 3));
    SET_VECTOR_ELT(list, 1, IntPtrToRVector(v2ptr, 4));
    SET_VECTOR_ELT(list, 2, CharPtrToRVector(v3ptr, 5));

    UNPROTECT(1);

    return list;
}

SEXP R_Rforce(
    SEXP splitFunctionIndex,
    SEXP interaction,
    SEXP designMatrixY,
    SEXP auxiliaryFeatures,
    SEXP varIDs,
    SEXP unitsOfCPIU,
    SEXP nTrees,
    SEXP maxDepth,
    SEXP minNodeSize,
    SEXP minGain,
    SEXP mtry,
    SEXP nsplits,
    SEXP seed)
{

    // convert R matrix to double**
    double **designMatrixY0 = RMatrixToDoublePtr(designMatrixY);
    double **auxiliaryFeatures0 = RMatrixToDoublePtr(auxiliaryFeatures);
    size_t nrow = Rf_nrows(designMatrixY);
    size_t ncolsDesign = Rf_ncols(designMatrixY);
    size_t ncolsAuxiliary = Rf_ncols(auxiliaryFeatures);

    // get variable IDs for categorical variables
    unsigned int *varIDs0 = (unsigned int *)INTEGER(varIDs);
    size_t nVars = Rf_length(varIDs);

    // get units of CPIU
    double *unitsOfCPIU0 = REAL(unitsOfCPIU);
    size_t nUnits = Rf_length(unitsOfCPIU);

    // get the hyperparameters for random forest
    size_t maxDepth0 = (size_t)INTEGER(maxDepth)[0];
    size_t minNodeSize0 = (size_t)INTEGER(minNodeSize)[0];
    double minGain0 = REAL(minGain)[0];
    size_t mtry0 = (size_t)INTEGER(mtry)[0];
    size_t nsplits0 = (size_t)INTEGER(nsplits)[0];
    size_t nTrees0 = (size_t)INTEGER(nTrees)[0];
    unsigned int seed0 = (unsigned int)INTEGER(seed)[0];

    // convert SEXP to C types, use case switch to select the function pointers later
    unsigned int splitFunctionIndex0 = (unsigned int)INTEGER(splitFunctionIndex)[0];

    fpSplitFunction splitFunction = NULL;
    fpLeafOutput leafOutputFunction = NULL;
    size_t lenOutput = 0;

    /* **********************
    initiate the phi array
    _treePhi is a 2D array of size nTrees * nUnits, a global variable split.h::34
    it is updated automatically in forest.c
    has to be allocated before calling the RandomForest function
    **************************/ 
    _treePhi = Allocate2DArray(nTrees0, nUnits);

    switch (splitFunctionIndex0)
    {
    case QuasiPoissonLikelihoodIndex:
        splitFunction = QuasiPoissonLikelihood;
        leafOutputFunction = LeafOutputInterval;
        lenOutput = nUnits;
        _pseudoRisk2 = 1L;
        _phi2 = 1L; 
        _longformat = 0L;
        _gee = 0L;
        break;
    case GEERuleIndex:
        splitFunction = GEERule;
        leafOutputFunction = LeafOutputInterval;
        lenOutput = nUnits;
        _pseudoRisk2 = 1L;
        _phi2 = 1L;
        _longformat = 0L;
        _gee = 1L;
        _padjust = "BH";
        _Rin = NULL;
        break;
    case PoissonLikelihoodIndex:
        splitFunction = PoissonLikelihood;
        leafOutputFunction = LeafOutput;
        lenOutput = 1L;
        _noPseudo = 1L;
        _longformat = 1L;
        _gee = 0L;
        break;
    default:
        Rf_error("Invalid splitFunctionIndex.");
    }

    // whether to account for interaction in gee split rule
    _interaction = (INTEGER(interaction)[0] == 1) ? 1: 0;
    ncolsDesign = ncolsDesign - lenOutput;

    // call the function
    RandomSurvivalForest *forest = RandomForest(
        splitFunction,
        leafOutputFunction,
        designMatrixY0,
        auxiliaryFeatures0,
        nrow,
        ncolsDesign,
        ncolsAuxiliary,
        varIDs0,
        nVars,
        unitsOfCPIU0,
        nUnits,
        lenOutput,
        maxDepth0,
        minNodeSize0,
        minGain0,
        mtry0,
        nsplits0,
        nTrees0,
        seed0);

    // rename the forest to forest0 for easy access, no need to free the memory as it is a shallow copy
    DecisionTreeNode **forest0 = forest->forest;

    // count the number of nodes in each tree
    size_t *nNodes = (size_t *)calloc(nTrees0, sizeof(size_t));
    size_t nNodesTotal = 0;
    for (size_t i = 0; i < nTrees0; i++)
    {
        nNodes[i] = CountNodes(forest0[i]);
        nNodesTotal += nNodes[i];
    }

    // forestMatrix: nNodesTotal * nDynamic(treeId, nodeId, status, size, leftDaughter, rightDaughter, splitVar, splitValue, prediction(lenOutput))
    double **forestMatrix = Allocate2DArray(nNodesTotal, 8 + lenOutput);

    size_t rowIndex = 0;
    for (size_t i = 0; i < nTrees0; i++)
    {
        SaveTreeArray(forest0[i], &rowIndex, &forestMatrix);
    }

    // return the forest to R
    // create an external pointer to the forest
    SEXP forestPtr = R_MakeExternalPtr(forest, R_NilValue, R_NilValue);
    R_RegisterCFinalizerEx(forestPtr, (R_CFinalizer_t)FreeSurvivalForest, TRUE);
    // bag matrix: nTrees * nrow
    SEXP bagMatrix = IntPtrToRMatrix(forest->bagMatrix, nTrees0, nrow);
    // tree_phi: nTrees * lenOutput
    SEXP treePhi = DoublePtrToRMatrix(_treePhi, nTrees0, nUnits);
    // convert the forestMatrix to R matrix
    SEXP forestMatrixR = DoublePtrToRMatrix(forestMatrix, nNodesTotal, 8 + lenOutput);
    // vimpStat: nVars * nTrees
    SEXP vimpStat = DoublePtrToRMatrix(forest->vimpPermuted, nTrees0, nVars);
    // predicted: nrow * lenOutput
    SEXP predicted = DoublePtrToRMatrix(forest->predicted, nrow, lenOutput);
    // oobPredicted: nrow * lenOutput
    SEXP oobPredicted = DoublePtrToRMatrix(forest->oobPredicted, nrow, lenOutput);
    // likelihoodsum: 1 * nTrees
    
    // free the memory
    // FreeSurvivalForest(forest);
    free(nNodes);
    Free2DArray(forestMatrix, nNodesTotal);

    // return the list
    SEXP list = PROTECT(Rf_allocVector(VECSXP, 7));
    SET_VECTOR_ELT(list, 0, bagMatrix);
    SET_VECTOR_ELT(list, 1, treePhi);
    SET_VECTOR_ELT(list, 2, forestMatrixR);
    SET_VECTOR_ELT(list, 3, vimpStat);
    SET_VECTOR_ELT(list, 4, predicted);
    SET_VECTOR_ELT(list, 5, oobPredicted);
    SET_VECTOR_ELT(list, 6, forestPtr); // Add as last element
    UNPROTECT(1);

    // set the names
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 7));
    SET_STRING_ELT(names, 0, Rf_mkChar("bagMatrix"));
    SET_STRING_ELT(names, 1, Rf_mkChar("treePhi"));
    SET_STRING_ELT(names, 2, Rf_mkChar("forestMatrix"));
    SET_STRING_ELT(names, 3, Rf_mkChar("vimpStat"));
    SET_STRING_ELT(names, 4, Rf_mkChar("predicted"));
    SET_STRING_ELT(names, 5, Rf_mkChar("oobPredicted"));
    SET_STRING_ELT(names, 6, Rf_mkChar("_external_forest_C_Ptr")); // Add name for the forest pointer
    Rf_setAttrib(list, R_NamesSymbol, names);
    UNPROTECT(1);

    return list;
}


SEXP R_ForestPredict(
    SEXP forestPtr,
    SEXP designMatrix)
{
    // convert the forest pointer to C pointer
    RandomSurvivalForest *forest = (RandomSurvivalForest *)R_ExternalPtrAddr(forestPtr);
    if (forest == NULL)
    {
        Rf_error("Invalid forest pointer.");
    }

    double **designMatrix0 = RMatrixToDoublePtr(designMatrix);
    size_t nrowsDesign = Rf_nrows(designMatrix);

    // call the ForestPredict function
    double **predicted = ForestPredict(
        forest->forest,
        forest->nTrees,
        designMatrix0,
        nrowsDesign,
        forest->lenOutput);

    // convert the predicted to R matrix
    SEXP predictedR = DoublePtrToRMatrix(predicted, nrowsDesign, forest->lenOutput);

    return predictedR;
}


SEXP R_SaveRforce(
    SEXP forestPtr,
    SEXP path)
{
    SaveSurvivalForest(
        (RandomSurvivalForest *)R_ExternalPtrAddr(forestPtr),
        CHAR(STRING_ELT(path, 0))
    );

    return R_NilValue;
}


SEXP R_LoadRforce(
    SEXP path)
{
    RandomSurvivalForest *forest = LoadSurvivalForest(CHAR(STRING_ELT(path, 0)));

    // create an external pointer to the forest
    SEXP forestPtr = R_MakeExternalPtr(forest, R_NilValue, R_NilValue);
    R_RegisterCFinalizerEx(forestPtr, (R_CFinalizer_t)FreeSurvivalForest, TRUE);

    return forestPtr;
}