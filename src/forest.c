#include "forest.h"

int _nPerms = 10;

RandomSurvivalForest *RandomForest(
    fpSplitFunction splitFunction,
    fpLeafOutput leafOutputFunction,
    double **designMatrixY,
    double **auxiliaryFeatures,
    size_t nrow,
    size_t ncolsDesign, // there are existing dummy variable coding for categorial variables, use varIDs to recover the vimpPermuted
    size_t ncolsAuxiliary,
    unsigned int *varIDs, // a vector of variable IDs, from 0 to (nVars -1). the consecutive repeated Ids indicates a single categorial variable when calculating vimpPermuted
    size_t nVars,
    double *unitsOfCPIU,
    size_t nUnits,
    size_t lenOutput,
    size_t maxDepth,
    size_t minNodeSize,
    double minGain,
    size_t mtry,
    size_t nsplits,
    size_t nTrees,
    unsigned int seed)
{
    DecisionTreeNode **forest = malloc(nTrees * sizeof(DecisionTreeNode *));

    double *patientIds = GetCol(auxiliaryFeatures, nrow, ncolsAuxiliary, 0);

    int **bootMatrix = BagMatrix(patientIds, nrow, nTrees, seed);

    double *likelihoodsum = (double *)calloc(nTrees, sizeof(double));
    double **vimp = Allocate2DArray(nTrees, ncolsDesign);

    free(patientIds);

    if (_longformat == 0)
    {
        // ignore the phi, phi = 1 all time
        if (_noPhi == 1 && _phi1 == 0 && _phi2 == 0 && _dynamicPhi == 0)
        {
            printf("Using phi = 1 for all splits in every tree\n\n");
            for (int i = 0; i < lenOutput; i++)
            {
                for (int j = 0; j < nTrees; j++)
                {
                    _treePhi[j][i] = 1.0;
                }
            }
        }

        // phi is calculated for the population
        if (_noPhi == 0 && _phi1 == 1 && _phi2 == 0 && _dynamicPhi == 0)
        {
            printf("Calculating phi once for the population\n\n");
            double *globalLambda = (double *)calloc(lenOutput, sizeof(double));
            double *globalYs = (double *)calloc(lenOutput, sizeof(double));
            double *globalRts = (double *)calloc(lenOutput, sizeof(double));
            double **globalYYs = Allocate2DArray(lenOutput, nrow);
            int *lenYY = (int *)calloc(lenOutput, sizeof(int));

            for (int i = 0; i < nrow; i++)
            {
                for (int j = 0; j < lenOutput; j++)
                {
                    if (fabs(auxiliaryFeatures[i][j + 1]) > 1e-9)
                    {
                        globalYs[j] += designMatrixY[i][ncolsDesign + j];
                        globalRts[j] += auxiliaryFeatures[i][j + 1];
                        globalYYs[j][lenYY[j]++] = designMatrixY[i][ncolsDesign + j] / sqrt(auxiliaryFeatures[i][j + 1]);
                    }
                }
            }

            double temp;

            for (int i = 0; i < lenOutput; i++)
            {
                globalLambda[i] = globalYs[i] / globalRts[i];
                temp = Var(globalYYs[i], lenYY[i]) / globalLambda[i];
                for (int j = 0; j < nTrees; j++)
                {
                    _treePhi[j][i] = temp;
                }
            }

            free(globalLambda);
            free(globalYs);
            free(globalRts);
            Free2DArray(globalYYs, lenOutput);
            free(lenYY);
        }
    }

    // phi is calculated for each tree
    if (_noPhi == 0 && _phi1 == 0 && _phi2 == 1 && _dynamicPhi == 0 && _longformat == 0)
    {
        printf("Calculating phi for each tree\n\n");
#pragma omp parallel for
        for (long i = 0; i < nTrees; i++)
        {
            DecisionTreeData *bootStrapData = BootStrapSample(designMatrixY,
                                                              auxiliaryFeatures,
                                                              nrow,
                                                              bootMatrix[i]);

            printf("Default: Calculating phi for each tree\n\n");
            double *globalLambda = (double *)calloc(lenOutput, sizeof(double));
            double *globalYs = (double *)calloc(lenOutput, sizeof(double));
            double *globalRts = (double *)calloc(lenOutput, sizeof(double));
            double **globalYYs = Allocate2DArray(lenOutput, nrow);
            int *lenYY = (int *)calloc(lenOutput, sizeof(int));

            for (int k = 0; k < nrow; k++)
            {
                for (int l = 0; l < lenOutput; l++)
                {
                    if (fabs(bootStrapData->auxiliaryFeatures[k][l + 1]) > 1e-9)
                    {
                        globalYs[l] += bootStrapData->designMatrixY[k][ncolsDesign + l];
                        globalRts[l] += bootStrapData->auxiliaryFeatures[k][l + 1];
                        globalYYs[l][lenYY[l]++] = bootStrapData->designMatrixY[k][ncolsDesign + l] / sqrt(bootStrapData->auxiliaryFeatures[k][l + 1]);
                    }
                }
            }

            for (int m = 0; m < lenOutput; m++)
            {
                globalLambda[m] = globalYs[m] / globalRts[m];
                _treePhi[i][m] = Var(globalYYs[m], lenYY[m]) / globalLambda[m];
            }

            free(globalLambda);
            free(globalYs);
            free(globalRts);
            Free2DArray(globalYYs, lenOutput);
            free(lenYY);

            free(bootStrapData[0].designMatrixY);
            free(bootStrapData[0].auxiliaryFeatures);
            free(bootStrapData[1].designMatrixY);
            free(bootStrapData[1].auxiliaryFeatures);
            free(bootStrapData);
        }
    }

    if (_verbose == 1)
    {
        printf("***************Phi****************\n");
        for (size_t nn = 0; nn < nTrees; nn++)
        {
            for (size_t mm = 0; mm < lenOutput; mm++)
            {
                printf("phi[%zu][%zu] is %f\n", nn, mm, _treePhi[nn][mm]);
            }
        }
    }

#pragma omp parallel for shared(_treePhi, _noPhi, _phi1, _phi2, _dynamicPhi, _noPseudo, _pseudoRisk1, _pseudoRisk2, _dynamicRisk, _k, _longformat, _gee, _verbose, _Rin, _padjust, _asympotic)
    for (long i = 0; i < nTrees; i++)
    {
        DecisionTreeData *bootStrapData = BootStrapSample(designMatrixY,
                                                          auxiliaryFeatures,
                                                          nrow,
                                                          bootMatrix[i]);

        double **designMatrixYIn = bootStrapData[0].designMatrixY;
        double **auxiliaryFeaturesIn = bootStrapData[0].auxiliaryFeatures;
        size_t nrowsIn = bootStrapData[0].nrows;

        double **designMatrixYOOB = bootStrapData[1].designMatrixY;
        double **auxiliaryFeaturesOOB = bootStrapData[1].auxiliaryFeatures;
        size_t nrowsOOB = bootStrapData[1].nrows;

        unsigned int treeSeed = seed + i;

        DecisionTreeNode *tree = Tree(
            splitFunction,
            leafOutputFunction,
            designMatrixYIn,
            auxiliaryFeaturesIn,
            nrowsIn,
            ncolsDesign,
            ncolsAuxiliary,
            unitsOfCPIU,
            nUnits,
            lenOutput,
            i,
            maxDepth,
            minNodeSize,
            minGain,
            mtry,
            nsplits,
            treeSeed);

        // calulated the likelihood for the original OOB data
        DecisionTreeNode *treeOOB = TreeOOB(
            tree,
            designMatrixYOOB,
            auxiliaryFeaturesOOB,
            nrowsOOB,
            ncolsDesign,
            ncolsAuxiliary,
            lenOutput);

        likelihoodsum[i] = TreeLikelihoodSum(treeOOB);

        // to store the cols to permute for each variable
        unsigned int *colsToPermute = (unsigned int *)calloc(nVars, sizeof(unsigned int));
        size_t lenColsToPermute = 0;
        // calculate the likelihood for the permutated OOB data at each column
        for (size_t k = 0; k < nVars; k++)
        {
            for (size_t kk = 0; kk < ncolsDesign; kk++)
            {
                if (varIDs[kk] == k)
                {
                    colsToPermute[lenColsToPermute++] = kk;
                }
            }

            double *tempLiklihoodPerm = (double *)calloc(_nPerms, sizeof(double));
            for (size_t l = 0; l < _nPerms; l++)
            {
                double **designMatrixYOOBPermuted = ColsPermute(
                    designMatrixYOOB,
                    nrowsOOB,
                    ncolsDesign + lenOutput,
                    colsToPermute,
                    lenColsToPermute,
                    treeSeed + l);

                DecisionTreeNode *treeOOBPermuted = TreeOOB(
                    tree,
                    designMatrixYOOBPermuted,
                    auxiliaryFeaturesOOB,
                    nrowsOOB,
                    ncolsDesign,
                    ncolsAuxiliary,
                    lenOutput);

                tempLiklihoodPerm[l] = TreeLikelihoodSum(treeOOBPermuted);
                Free2DArray(designMatrixYOOBPermuted, nrowsOOB);
                FreeTree(treeOOBPermuted);
            }

            vimp[i][k] = likelihoodsum[i] - Mean(tempLiklihoodPerm, _nPerms);
            free(tempLiklihoodPerm);
            lenColsToPermute = 0;
        }
        free(colsToPermute);

#pragma omp critical
        {
            forest[i] = tree;
        }

        free(designMatrixYIn);
        free(auxiliaryFeaturesIn);
        free(designMatrixYOOB);
        free(auxiliaryFeaturesOOB);
        free(bootStrapData);
        FreeTree(treeOOB);
    }

    double **oobPredicted = InternalForestOOBPredict(
        forest,
        nTrees,
        designMatrixY,
        nrow,
        lenOutput,
        bootMatrix);

    double **predicted = ForestPredict(
        forest,
        nTrees,
        designMatrixY,
        nrow,
        lenOutput);

    double *vimpStat = (double *)calloc(ncolsDesign, sizeof(double));
    double *vimpFreq = (double *)calloc(ncolsDesign, sizeof(double));

    VimpForest(
        forest,
        nTrees,
        &vimpStat,
        &vimpFreq);

    // calculate the mse, both inbag and oob
    double **oobPredictedNEvents = Allocate2DArray(nrow, lenOutput);
    double **inbagPredictedNEvents = Allocate2DArray(nrow, lenOutput);

    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < lenOutput; j++)
        {
            // predicted lambda * (pseudo) risk time
            oobPredictedNEvents[i][j] = oobPredicted[i][j] * auxiliaryFeatures[i][j + 1];
            inbagPredictedNEvents[i][j] = predicted[i][j] * auxiliaryFeatures[i][j + 1];
        }
    }

    double oobMSE = 0;
    double inbagMSE = 0;

    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < lenOutput; j++)
        {
            oobMSE += pow(designMatrixY[i][ncolsDesign + j] - oobPredictedNEvents[i][j], 2) / (nrow * lenOutput);
            inbagMSE += pow(designMatrixY[i][ncolsDesign + j] - inbagPredictedNEvents[i][j], 2) / (nrow * lenOutput);
        }
    }

    Free2DArray(oobPredictedNEvents, nrow);
    Free2DArray(inbagPredictedNEvents, nrow);

    printf("\n\n##############################\n\n");
    printf("VimpFreq: \n");
    PrintArrayDouble(vimpFreq, ncolsDesign);
    printf("VimpStat: \n");
    PrintArrayDouble(vimpStat, ncolsDesign);
    printf("oobMSE: %.5f, inbagMSE: %0.5f\n", oobMSE, inbagMSE);

    RandomSurvivalForest *out = malloc(sizeof(RandomSurvivalForest));

    // copy of unitsOfCPIU, so when free the forest, the unitsOfCPIU will not be freed but the copy
    double *unitsOfCPIUCopy = (double *)calloc(nUnits, sizeof(double));
    memcpy(unitsOfCPIUCopy, unitsOfCPIU, nUnits * sizeof(double));

    out[0] = (RandomSurvivalForest){
        nrow,
        ncolsDesign,
        varIDs,
        nVars,
        nTrees,
        unitsOfCPIUCopy,
        nUnits,
        lenOutput,
        maxDepth,
        minNodeSize,
        minGain,
        mtry,
        nsplits,
        seed,

        predicted,
        oobPredicted,

        oobMSE,
        inbagMSE,

        vimpStat,
        vimpFreq,
        likelihoodsum,
        vimp,

        bootMatrix,

        forest};

    // printf("seed is %d\n", out[0].seed);

    return out;
}

double **InternalForestOOBPredict(
    DecisionTreeNode **forest,
    size_t nTrees,
    double **designMatrix,
    size_t nrowsDesign,
    size_t lenOutput,
    int **BagMatrix)
{
    double ***prediction = Allocate3DArray(nrowsDesign, nTrees, lenOutput);

#pragma omp parallel
    {
#pragma omp for nowait
        for (int i = 0; i < nrowsDesign; i++)
        {
            for (int j = 0; j < nTrees; j++)
            {
                if (BagMatrix[j][i] != 0)
                {
                    double *temp = TreePredict(
                        forest[j],
                        designMatrix[i]);
                    memcpy(prediction[i][j], temp, lenOutput * sizeof(double));
                }
            }
        }
    }

    // get total # of oob trees for each row of design matrix
    int *nOOBTrees = (int *)calloc(nrowsDesign, sizeof(int));

    for (int i = 0; i < nrowsDesign; i++)
    {
        for (int j = 0; j < nTrees; j++)
        {
            if (BagMatrix[j][i] != 0)
            {
                nOOBTrees[i] += 1;
            }
        }
    }

    double **out = Allocate2DArray(nrowsDesign, lenOutput);

    for (int i = 0; i < nrowsDesign; i++)
    {
        for (int j = 0; j < nTrees; j++)
        {
            for (int k = 0; k < lenOutput; k++)
            {
                if (BagMatrix[j][i] != 0)
                {
                    out[i][k] += prediction[i][j][k] / nOOBTrees[i];
                    // printf("prediction[%d][%d][%d] is %f\n", i, j, k, prediction[i][j][k]);
                }
            }
        }
    }

    Free3DArray(prediction, nrowsDesign, nTrees);
    free(nOOBTrees);

    return out;
}

double **ForestPredict(
    DecisionTreeNode **forest,
    size_t nTrees,
    double **designMatrix,
    size_t nrowsDesign,
    size_t lenOutput)
{
    double ***prediction = Allocate3DArray(nrowsDesign, nTrees, lenOutput);

#pragma omp parallel
    {
#pragma omp for nowait
        for (int i = 0; i < nrowsDesign; i++)
        {
            for (int j = 0; j < nTrees; j++)
            {
                double *temp = TreePredict(
                    forest[j],
                    designMatrix[i]);
                memcpy(prediction[i][j], temp, lenOutput * sizeof(double));
            }
        }
    }

    double **out = Allocate2DArray(nrowsDesign, lenOutput);

    for (int i = 0; i < nrowsDesign; i++)
    {
        for (int j = 0; j < nTrees; j++)
        {
            for (int k = 0; k < lenOutput; k++)
            {
                out[i][k] += prediction[i][j][k] / nTrees;
            }
        }
    }

    Free3DArray(prediction, nrowsDesign, nTrees);
    return out;
}

double **SurvivalForestPredict(
    RandomSurvivalForest *survivalForest,
    double **designMatrix,
    size_t nrowsDesign)
{
    double **predicted = ForestPredict(
        survivalForest->forest,
        survivalForest->nTrees,
        designMatrix,
        nrowsDesign,
        survivalForest->lenOutput);

    return predicted;
}

void VimpForest(
    DecisionTreeNode **forest,
    size_t nTrees,
    double **vimpStat,
    double **vimpFreq)
{
    for (int i = 0; i < nTrees; i++)
    {
        VimpTree(forest[i], vimpStat, vimpFreq);
    }
}

void FreeSurvivalForest(RandomSurvivalForest *forest)
{
    for (int i = 0; i < forest->nTrees; i++)
    {
        FreeTree(forest->forest[i]);
    }
    free(forest->forest);
    Free2DArray(forest->predicted, forest->nrowsDesign);
    Free2DArray(forest->oobPredicted, forest->nrowsDesign);
    free(forest->vimpStat);
    free(forest->vimpFreq);
    free(forest->likelihoodsum);
    Free2DArrayInt(forest->bagMatrix, forest->nTrees);
    Free2DArray(forest->vimpPermuted, forest->nTrees);
    // as unitsOfCPIU is a pointer from outside, it should not be freed here, make it easy to debug
    if (forest->unitsOfCPIU != NULL)
        free(forest->unitsOfCPIU);
    free(forest);
}

void SaveForest(DecisionTreeNode **forest,
                size_t nTrees,
                size_t lenOutput,
                char *path)
{
    char outputFolderTree[300];
    char outputFolderDot[300];
    sprintf(outputFolderTree, "%s/%s", path, "trees");
    sprintf(outputFolderDot, "%s/%s", path, "dots");

    if (access(outputFolderTree, F_OK) == -1)
    {
        if (mkdir(outputFolderTree, 0777) == -1)
        {
            printf("Could not create folder %s!\n", outputFolderTree);
            exit(1);
        }
    }

    if (access(outputFolderDot, F_OK) == -1)
    {
        if (mkdir(outputFolderDot, 0777) == -1)
        {
            printf("Could not create folder %s!\n", outputFolderDot);
            exit(1);
        }
    }
    char *pathTree = malloc(300 * sizeof(char));
    char *pathDot = malloc(300 * sizeof(char));

    for (int i = 0; i < nTrees; i++)
    {
        sprintf(pathTree, "%s/tree_%d.tdb", outputFolderTree, i);
        FILE *file = fopen(pathTree, "wb");
        SaveTree(forest[i], file);
        fclose(file);

        sprintf(pathDot, "%s/tree_%d.dot", outputFolderDot, i);
        WriteTreeDotFile(forest[i], pathDot, 0);
    }

    free(pathTree);
    free(pathDot);
}

void SaveSurvivalForest(RandomSurvivalForest *forest,
                        char *path)
{
    // save forest
    printf("seed is %ld\n", forest->seed);
    char *pathForest = malloc(1024 * sizeof(char));
    sprintf(pathForest, "%s/%s_tree%zu_depth%zu_mtry%zu_Gain%.3f_nsplits%zu_nodesize%zu_seed%d",
            path,
            "forest",
            forest->nTrees,
            forest->maxDepth,
            forest->mtry,
            forest->minGain,
            forest->nsplits,
            forest->minNodeSize,
            forest->seed);

    MkdirRecursive(pathForest);

    printf("Saving forest to %s\n", pathForest);

    SaveForest(forest->forest, forest->nTrees, forest->lenOutput, pathForest);

    // save forest parameters
    char *pathParameters = malloc(300 * sizeof(char));
    sprintf(pathParameters, "%s/%s", pathForest, "parameters.tdb");
    FILE *fileParameters = fopen(pathParameters, "wb");

    fprintf(fileParameters, "%ld,%ld,%ld,%ld\n", forest->nrowsDesign, forest->ncolsDesign, forest->nVars, forest->nTrees);

    if (forest->unitsOfCPIU != NULL)
    {
        // 1 indicating that wide format is used
        fprintf(fileParameters, "%d,%ld,", 1, forest->nUnits);
        for (int i = 0; i < (forest->nUnits - 1); i++)
        {
            fprintf(fileParameters, "%.5f,", forest->unitsOfCPIU[i]);
        }
        fprintf(fileParameters, "%.5f\n", forest->unitsOfCPIU[forest->nUnits - 1]);
    }
    else
    {
        fprintf(fileParameters, "%d,%ld,\n", 0, forest->nUnits);
    }

    fprintf(
        fileParameters,
        "%ld,%ld,%ld,%.5f,%ld,%ld,%d\n",
        forest->lenOutput,
        forest->maxDepth,
        forest->minNodeSize,
        forest->minGain,
        forest->mtry,
        forest->nsplits,
        forest->seed);

    fprintf(fileParameters,
            "%.5f, %.5f\n",
            forest->oobMSE,
            forest->inbagMSE);

    fclose(fileParameters);
    free(pathParameters);

    // save predicted results of training data
    char *pathPredicted = malloc(300 * sizeof(char));
    char *pathPredictedOob = malloc(300 * sizeof(char));
    char *pathLikelihoodSum = malloc(300 * sizeof(char));
    char *pathvimpPermuted = malloc(300 * sizeof(char));

    sprintf(pathPredicted, "%s/%s", pathForest, "predicted.tdb");
    sprintf(pathPredictedOob, "%s/%s", pathForest, "predicted_oob.tdb");
    sprintf(pathLikelihoodSum, "%s/%s", pathForest, "likelihoodSum.tdb");
    sprintf(pathvimpPermuted, "%s/%s", pathForest, "vimpPermuted.tdb");

    FILE *filePredicted = fopen(pathPredicted, "wb");
    FILE *filePredictedOob = fopen(pathPredictedOob, "wb");
    FILE *fileLikelihoodSum = fopen(pathLikelihoodSum, "wb");
    FILE *filevimpPermuted = fopen(pathvimpPermuted, "wb");

    WriteCSV(forest->predicted, filePredicted, forest->nrowsDesign, forest->lenOutput);
    WriteCSV(forest->oobPredicted, filePredictedOob, forest->nrowsDesign, forest->lenOutput);
    WriteCSV(&(forest->likelihoodsum), fileLikelihoodSum, 1, forest->nTrees);
    WriteCSV(forest->vimpPermuted, filevimpPermuted, forest->nTrees, forest->nVars);

    fclose(filePredicted);
    fclose(filePredictedOob);
    fclose(fileLikelihoodSum);
    fclose(filevimpPermuted);

    free(pathPredicted);
    free(pathPredictedOob);
    free(pathLikelihoodSum);
    free(pathvimpPermuted);

    // save vimp
    char *pathVimp = malloc(300 * sizeof(char));
    sprintf(pathVimp, "%s/%s", pathForest, "vimp.tdb");
    FILE *fileVimp = fopen(pathVimp, "wb");

    WriteCSV(&forest->vimpStat, fileVimp, 1, forest->ncolsDesign);
    WriteCSV(&forest->vimpFreq, fileVimp, 1, forest->ncolsDesign);

    fclose(fileVimp);
    free(pathVimp);

    // save bag matrix
    char *pathBagMatrix = malloc(300 * sizeof(char));
    sprintf(pathBagMatrix, "%s/%s", pathForest, "bagMatrix.tdb");
    FILE *fileBagMatrix = fopen(pathBagMatrix, "wb");

    WriteCSVInt(forest->bagMatrix, fileBagMatrix, forest->nTrees, forest->nrowsDesign);

    fclose(fileBagMatrix);
    free(pathBagMatrix);
    free(pathForest);
}

DecisionTreeNode **LoadForest(char *path,
                              size_t nTrees)
{
    if (access(path, F_OK) == -1)
    {
        printf("Input folder %s does not exist!\n", path);
        return 0;
    }

    DecisionTreeNode **forest = malloc(nTrees * sizeof(DecisionTreeNode *));

    char *pathTree = malloc(300 * sizeof(char));

    for (int i = 0; i < nTrees; i++)
    {
        sprintf(pathTree, "%s/trees/tree_%d.tdb", path, i);
        FILE *file = fopen(pathTree, "rb");
        forest[i] = LoadTree(file);
        fclose(file);
    }

    free(pathTree);

    return forest;
}

RandomSurvivalForest *LoadSurvivalForest(char *path)
{
    RandomSurvivalForest *forest = (RandomSurvivalForest *)malloc(sizeof(RandomSurvivalForest));

    // load forest parameters
    char *pathParameters = malloc(300 * sizeof(char));

    sprintf(pathParameters, "%s/%s", path, "parameters.tdb");
    FILE *fileParameters = fopen(pathParameters, "rb");

    if (fileParameters == NULL)
    {
        printf("Cannot open file %s\n", pathParameters);
        return NULL;
    }

    printf("Loading forest parameters from %s\n", pathParameters);
    fscanf(fileParameters, "%ld,%ld,%ld,%ld\n", &forest->nrowsDesign, &forest->ncolsDesign, &forest->nVars, &forest->nTrees);
    printf("nrowsDesign: %ld, ncolsDesign: %ld, nVars: %ld, nTrees: %ld\n", forest->nrowsDesign, forest->ncolsDesign, forest->nVars, forest->nTrees);

    int hasUnitsOfCPIU;
    fscanf(fileParameters, "%d,%ld,", &hasUnitsOfCPIU, &forest->nUnits);
    if (hasUnitsOfCPIU)
    {
        double *unitsOfCPIU = malloc(forest->nrowsDesign * sizeof(double));

        for (int i = 0; i < (forest->nUnits - 1); i++)
        {
            fscanf(fileParameters, "%lf,", &unitsOfCPIU[i]);
        }
        fscanf(fileParameters, "%lf\n", &unitsOfCPIU[forest->nUnits - 1]);
        forest->unitsOfCPIU = unitsOfCPIU;
    }
    else
    {
        fscanf(fileParameters, "\n"); // skip the rest of the line
        forest->unitsOfCPIU = NULL;
    }

    fscanf(
        fileParameters,
        "%ld,%ld,%ld,%lf,%ld,%ld,%d\n",
        &forest->lenOutput,
        &forest->maxDepth,
        &forest->minNodeSize,
        &forest->minGain,
        &forest->mtry,
        &forest->nsplits,
        &forest->seed);

    fscanf(fileParameters,
           "%lf, %lf\n",
           &forest->oobMSE,
           &forest->inbagMSE);

    // printf("lenOutput: %ld, maxDepth: %ld, minNodeSize: %ld, minGain: %lf, mtry: %ld, nsplits: %ld, seed: %d\n", forest->lenOutput, forest->maxDepth, forest->minNodeSize, forest->minGain, forest->mtry, forest->nsplits, forest->seed);
    fclose(fileParameters);
    free(pathParameters);

    // load predicted results of training data
    char *pathPredicted = malloc(300 * sizeof(char));
    char *pathPredictedOob = malloc(300 * sizeof(char));
    char *pathLikelihoodSum = malloc(300 * sizeof(char));
    char *pathvimpPermuted = malloc(300 * sizeof(char));

    sprintf(pathPredicted, "%s/%s", path, "predicted.tdb");
    sprintf(pathPredictedOob, "%s/%s", path, "predicted_oob.tdb");
    sprintf(pathLikelihoodSum, "%s/%s", path, "likelihoodSum.tdb");
    sprintf(pathvimpPermuted, "%s/%s", path, "vimpPermuted.tdb");

    double **predicted = Allocate2DArray(forest->nrowsDesign, forest->lenOutput);
    double **predictedOob = Allocate2DArray(forest->nrowsDesign, forest->lenOutput);
    double **likelihoodSum = Allocate2DArray(1, forest->nTrees);
    double **vimpPermuted = Allocate2DArray(forest->nTrees, forest->nVars);

    // printf("nrowsDesign: %ld, nVars: %ld, nTrees: %ld\n", forest->nrowsDesign, forest->nVars, forest->nTrees);
    ReadCSV(pathPredicted, &predicted, 0);
    ReadCSV(pathPredictedOob, &predictedOob, 0);
    ReadCSV(pathLikelihoodSum, &likelihoodSum, 0);
    ReadCSV(pathvimpPermuted, &vimpPermuted, 0);

    forest->predicted = predicted;
    forest->oobPredicted = predictedOob;
    forest->likelihoodsum = (double *)malloc(forest->nTrees * sizeof(double));
    memcpy(forest->likelihoodsum, likelihoodSum[0], forest->nTrees * sizeof(double));
    forest->vimpPermuted = vimpPermuted;

    free(pathPredicted);
    free(pathPredictedOob);
    free(pathLikelihoodSum);
    free(pathvimpPermuted);

    Free2DArray(likelihoodSum, 1);

    // load vimp
    char *pathVimp = malloc(300 * sizeof(char));
    sprintf(pathVimp, "%s/%s", path, "vimp.tdb");

    double **vimp = Allocate2DArray(2, forest->ncolsDesign);
    ReadCSV(pathVimp, &vimp, 0);

    forest->vimpStat = vimp[0];
    forest->vimpFreq = vimp[1];
    free(pathVimp);
    free(vimp);

    // load bag matrix
    char *pathBagMatrix = malloc(300 * sizeof(char));
    sprintf(pathBagMatrix, "%s/%s", path, "bagMatrix.tdb");

    int **bagMatrix = Allocate2DArrayInt(forest->nTrees, forest->nrowsDesign);
    // printf("nTrees: %ld, nrowsDesign: %ld\n", forest->nTrees, forest->nrowsDesign);
    ReadCSVInt(pathBagMatrix, &bagMatrix, 0);

    forest->bagMatrix = bagMatrix;
    free(pathBagMatrix);

    // load forest

    forest->forest = LoadForest(path, forest->nTrees);

    return forest;
}

void PrintForestPaths(RandomSurvivalForest *forest, FILE *file)
{
    if (forest == NULL)
    {
        printf("Forest is NULL!\n");
        return;
    }

    for (int i = 0; i < forest->nTrees; i++)
    {
        PrintTreeNodes(forest->forest[i], file);
        fprintf(file, "\n\n");
    }
}
