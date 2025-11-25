#include "tree.h"

DecisionTreeNode *EmptyNode(long treeId, long *nodeId)
{
    DecisionTreeNode *node = (DecisionTreeNode *)malloc(sizeof(DecisionTreeNode));

    if (node == NULL)
    {
        PRINT_LOCATION();
        printf("Error: memory allocation failed in EmptyNode()!");
        exit(1);
    }

    node->nodeId = *nodeId;
    node->treeId = treeId;
    node->flag = 0;
    node->splitIndex = -1;
    node->splitValue = 0.0;
    node->splitStat = NULL;
    node->sizeLR = NULL;
    node->output = NULL;
    node->lenOutput = 0;

    node->leftChild = NULL;
    node->rightChild = NULL;

    (*nodeId)++;

    return node;
}

void GrowTree(
    fpSplitFunction splitFunction,
    fpLeafOutput leafOutputFunction,
    DecisionTreeNode *tree,
    double **designMatrixY,     // only a copy of pointer to rows, not entire data matrix
    double **auxiliaryFeatures, // only a copy of pointer to rows, not entire data matrix
    size_t nrows,               // nrow for design matrix
    size_t ncolsDesign,         // ncol for design matrix, including Y
    size_t ncolsAuxiliary,      // ncol for auxiliaryFeature matrix
    double *unitsOfCPIU,        // length of each CPIU
    size_t nUnits,              // number of CPIU
    size_t lenOutput,           // length of output
    long treeId,                // current tree id
    long *nodeId,               // current node id
    size_t depth,               // current depth
    size_t maxDepth,            // max depth
    size_t minNodeSize,         // min size of leaf node
    double minGain,             // min gain of split
    size_t mtry,                // mtry variables when splitting
    size_t nsplits,             // maximum attempts of each variable when splitting
    unsigned int seed)
{
    // if the node is a leaf node, return
    if (tree->flag == 1)
    {
        return;
    }

    // check if depth is larger than maxDepth
    if (depth > maxDepth)
    {
        tree->flag = 1;
        tree->splitIndex = -1;
        tree->output = leafOutputFunction(
            nrows,
            ncolsDesign,
            ncolsAuxiliary,
            designMatrixY,
            auxiliaryFeatures,
            unitsOfCPIU,
            treeId,
            nUnits,
            lenOutput);
        tree->lenOutput = lenOutput;
        free(designMatrixY);
        free(auxiliaryFeatures);
        size_t *nLR = (size_t *)calloc(3, sizeof(size_t));
        nLR[2] = nrows;

        tree->sizeLR = nLR;
        return;
    }

    // find the best split
    SplitPoints *bestSplit = FindBestSplit(
        splitFunction,
        leafOutputFunction,
        designMatrixY,
        auxiliaryFeatures,
        unitsOfCPIU,
        treeId,
        nUnits,
        lenOutput,
        mtry,
        nsplits,
        nrows,
        ncolsDesign,
        ncolsAuxiliary,
        minNodeSize,
        seed);

    // if the split is not found or not good enough, make the node a leaf node
    if (bestSplit->splitIndex == -1 ||
        bestSplit->splitStat[0] < minGain)
    {
        tree->flag = 1;
        tree->splitIndex = -1;
        tree->output = leafOutputFunction(
            nrows,
            ncolsDesign,
            ncolsAuxiliary,
            designMatrixY,
            auxiliaryFeatures,
            unitsOfCPIU,
            treeId,
            nUnits,
            lenOutput);
        tree->lenOutput = lenOutput;
        free(designMatrixY);
        free(auxiliaryFeatures);
        free(bestSplit->splitStat);
        free(bestSplit);
        size_t *nLR = (size_t *)calloc(3, sizeof(size_t));
        nLR[2] = nrows;

        tree->sizeLR = nLR;
        return;
    }
    else
    {
        // populate the best split to tree node
        tree->splitIndex = bestSplit->splitIndex;
        tree->splitValue = bestSplit->splitValue;
        tree->splitStat = bestSplit->splitStat;

        // free the memory
        free(bestSplit);

        // split the data based on the best split
        DecisionTreeData *dataSplits = SplitDataset(
            tree->splitIndex,
            tree->splitValue,
            nrows,
            designMatrixY,
            auxiliaryFeatures);

        // rename the data splits
        double **designMatrixYLeft = dataSplits[0].designMatrixY;
        double **designMatrixYRight = dataSplits[1].designMatrixY;
        double **auxiliaryFeaturesLeft = dataSplits[0].auxiliaryFeatures;
        double **auxiliaryFeaturesRight = dataSplits[1].auxiliaryFeatures;
        size_t nLeft = dataSplits[0].nrows;
        size_t nRight = dataSplits[1].nrows;

        size_t *nLR = (size_t *)calloc(3, sizeof(size_t));
        nLR[0] = nLeft;
        nLR[1] = nRight;
        nLR[2] = nrows;

        tree->sizeLR = nLR;

        // free the memory
        free(dataSplits);

        // check if daughter nodes are too small
        if (nLeft < minNodeSize || nRight < minNodeSize)
        {
            tree->flag = 1;
            tree->splitIndex = -1;
            tree->output = leafOutputFunction(
                nrows,
                ncolsDesign,
                ncolsAuxiliary,
                designMatrixY,
                auxiliaryFeatures,
                unitsOfCPIU,
                treeId,
                nUnits,
                lenOutput);
            tree->lenOutput = lenOutput;
            free(designMatrixYLeft);
            free(designMatrixYRight);
            free(auxiliaryFeaturesLeft);
            free(auxiliaryFeaturesRight);
            free(designMatrixY);
            free(auxiliaryFeatures);
            return;
        }
        else
        {
            // free the memory
            free(designMatrixY);
            free(auxiliaryFeatures);

            // create left child node
            tree->leftChild = EmptyNode(treeId, nodeId);
            GrowTree(
                splitFunction,
                leafOutputFunction,
                tree->leftChild,
                designMatrixYLeft,
                auxiliaryFeaturesLeft,
                nLeft,
                ncolsDesign,
                ncolsAuxiliary,
                unitsOfCPIU,
                nUnits,
                lenOutput,
                treeId,
                nodeId,
                depth + 1,
                maxDepth,
                minNodeSize,
                minGain,
                mtry,
                nsplits,
                tree->sizeLR[0]);

            // create right child node
            tree->rightChild = EmptyNode(treeId, nodeId);
            GrowTree(
                splitFunction,
                leafOutputFunction,
                tree->rightChild,
                designMatrixYRight,
                auxiliaryFeaturesRight,
                nRight,
                ncolsDesign,
                ncolsAuxiliary,
                unitsOfCPIU,
                nUnits,
                lenOutput,
                treeId,
                nodeId,
                depth + 1,
                maxDepth,
                minNodeSize,
                minGain,
                mtry,
                nsplits,
                tree->sizeLR[1]);
        }
    }
}

DecisionTreeNode *Tree(
    fpSplitFunction splitFunction,
    fpLeafOutput leafOutputFunction,
    double **designMatrixY,
    double **auxiliaryFeatures,
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double *unitsOfCPIUs,
    size_t nUnits,
    size_t lenOutput,
    long treeId,
    size_t maxDepth,
    size_t minNodeSize,
    double minGain,
    size_t mtry,
    size_t nsplits,
    unsigned int seed)
{

    // create the root node
    long nodeId = 0;
    DecisionTreeNode *root = EmptyNode(treeId, &nodeId);

    // copy the design matrix first
    double **designMatrixYCopy = (double **)malloc(nrows * sizeof(double *));
    double **auxiliaryFeaturesCopy = (double **)malloc(nrows * sizeof(double *));

    if (_noPseudo == 1 || _pseudoRisk1 == 1 || _longformat == 1)
    {
        // printf("Keep the input risk time for tree %ld...\n", treeId);
        /****************************Use the global pseudo-risk time as-is****************/
        // no pseudo risk time and population level pseudo risk time are handled by the user input
        // pseudo risk time at split is handled in split.c
        // pseudo risk time at tree-level is handled in current else condition
        for (int i = 0; i < nrows; i++)
        {
            designMatrixYCopy[i] = designMatrixY[i];
            auxiliaryFeaturesCopy[i] = auxiliaryFeatures[i];
        }

        GrowTree(
            splitFunction,
            leafOutputFunction,
            root,
            designMatrixYCopy,
            auxiliaryFeaturesCopy,
            nrows,
            ncolsDesign,
            ncolsAuxiliary,
            unitsOfCPIUs,
            nUnits,
            lenOutput,
            treeId,
            &nodeId,
            0,
            maxDepth,
            minNodeSize,
            minGain,
            mtry,
            nsplits,
            seed);

        return root;
    }
    else
    {
        /**************recalculate pseudo-risk time for each tree*********************/
        // printf("Default: recalculate pseudo-risk time for tree %ld...\n", treeId);
        // no pseudo risk time and population level pseudo risk time are handled by the user input
        // pseudo risk time at split is handled in split.c
        // pseudo risk time at tree-level is handled in current else condition
        double **pseudoRiskTimes = Allocate2DArray(nrows, nUnits + 3);
        double *X = GetCol(auxiliaryFeatures, nrows, ncolsAuxiliary, ncolsAuxiliary - 2);
        double *status = GetCol(auxiliaryFeatures, nrows, ncolsAuxiliary, ncolsAuxiliary - 1);
        ElementStruct *uniqueStatus = Unique(status, nrows);
        if (uniqueStatus->nElements > 2)
        {
            PRINT_LOCATION();
            printf("Error: status should only have at most two unique values.\n The current unique values are:\n");
            PrintArrayDouble(uniqueStatus->Elements, uniqueStatus->nElements);
            exit(1);
        }
        free(uniqueStatus->Elements);
        free(uniqueStatus);
        double *statusInverse = (double *)malloc(nrows * sizeof(double));

        for (int i = 0; i < nrows; i++)
        {
            if (fabs(status[i]) < 1e-6)
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

        KMResult *wtPatient; // placeholder for inverse probability of censoring weight w.r.t time of interest

        for (int i = 0; i < nrows; i++)
        {
            wtPatient = Wt(Gt, X[i], status[i]);
            pseudoRiskTimes[i][0] = auxiliaryFeatures[i][0]; // patient id
            pseudoRiskTimes[i][nUnits + 1] = X[i];           // time of interest
            pseudoRiskTimes[i][nUnits + 2] = status[i];      // censoring status
            for (int j = 0; j < nUnits; j++)
            {
                pseudoRiskTimes[i][j + 1] = PseudoRiskTime(wtPatient, breakpoints[j], breakpoints[j + 1]);
            }
            FreeKMResult(wtPatient);
        }
        FreeKMResult(Gt);
        free(X);
        free(status);
        for (int i = 0; i < nrows; i++)
        {
            designMatrixYCopy[i] = designMatrixY[i];
            auxiliaryFeaturesCopy[i] = pseudoRiskTimes[i];
        }

        GrowTree(
            splitFunction,
            leafOutputFunction,
            root,
            designMatrixYCopy,
            auxiliaryFeaturesCopy,
            nrows,
            ncolsDesign,
            ncolsAuxiliary,
            unitsOfCPIUs,
            nUnits,
            lenOutput,
            treeId,
            &nodeId,
            0,
            maxDepth,
            minNodeSize,
            minGain,
            mtry,
            nsplits,
            seed);

        Free2DArray(pseudoRiskTimes, nrows);
        free(breakpoints);
        return root;
    }
}

DecisionTreeNode *TreeCopy(DecisionTreeNode *root)
{
    if (root == NULL)
    {
        return NULL;
    }

    DecisionTreeNode *copy = (DecisionTreeNode *)malloc(sizeof(DecisionTreeNode));
    copy->treeId = root->treeId;
    copy->nodeId = root->nodeId;
    copy->flag = root->flag;
    copy->splitIndex = root->splitIndex;
    copy->splitValue = root->splitValue;
    copy->lenOutput = root->lenOutput;
    if (root->flag == 0)
    {
        copy->sizeLR = (size_t *)malloc(2 * sizeof(size_t));
        copy->sizeLR[0] = root->sizeLR[0];
        copy->sizeLR[1] = root->sizeLR[1];
        copy->splitStat = (double *)malloc(4 * sizeof(double));
        copy->splitStat[0] = root->splitStat[0];
        copy->splitStat[1] = root->splitStat[1];
        copy->splitStat[2] = root->splitStat[2];
        copy->splitStat[3] = root->splitStat[3];
        copy->output = NULL;
    }
    else
    {
        copy->sizeLR = NULL;
        copy->splitStat = NULL;
        copy->output = (double *)malloc(copy->lenOutput * sizeof(double));
        for (size_t i = 0; i < copy->lenOutput; i++)
        {
            copy->output[i] = root->output[i];
        }
    }

    copy->leftChild = TreeCopy(root->leftChild);
    copy->rightChild = TreeCopy(root->rightChild);

    return copy;
}

void TreeOOBInternal(
    DecisionTreeNode *rootCopy,
    double **designMatrixY,     // the out-of-bag design matrix
    double **auxiliaryFeatures, // the out-of-bag auxiliary features, pseudo-risk estimation should be done in forest step before feeding into this function, similar with GrowTree()
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    size_t lenOutput)
{
    if (rootCopy->flag == 1)
    {
        double temp = 0.0f, temp1 = 0.0f;
        for (size_t i = 0; i < nrows; i++)
        {
            for (size_t j = 0; j < lenOutput; j++)
            {
                if (fabs(auxiliaryFeatures[i][j + 1]) < 1e-9 || fabs((rootCopy->output)[j]) < 1e-9)
                    continue;
                if (_verbose == 2)
                {
                    printf("designMatrixY[%zu][%zu] = %lf; ", i, j, designMatrixY[i][j]);
                    printf("auxiliaryFeatures[%zu][%zu] = %lf; ", i, j + 1, auxiliaryFeatures[i][j + 1]);
                    printf("rootCopy->output[%zu] = %lf; ", j, (rootCopy->output)[j]);
                    printf("single likelihood = %lf\n", temp1);
                }
                temp1 = 1.0f / _treePhi[rootCopy->treeId][j] * (designMatrixY[i][ncolsDesign + j] * log((rootCopy->output)[j] * auxiliaryFeatures[i][j + 1]) - (rootCopy->output)[j] * auxiliaryFeatures[i][j + 1]);
                temp += temp1;
            }
        }
        rootCopy->splitStat = (double *)malloc(1 * sizeof(double));
        rootCopy->splitStat[0] = temp;
        return;
    }

    DecisionTreeData *dataSplits = SplitDataset(
        rootCopy->splitIndex,
        rootCopy->splitValue,
        nrows,
        designMatrixY,
        auxiliaryFeatures);

    // rename the data splits
    double **designMatrixYLeft = dataSplits[0].designMatrixY;
    double **designMatrixYRight = dataSplits[1].designMatrixY;
    double **auxiliaryFeaturesLeft = dataSplits[0].auxiliaryFeatures;
    double **auxiliaryFeaturesRight = dataSplits[1].auxiliaryFeatures;
    size_t nLeft = dataSplits[0].nrows;
    size_t nRight = dataSplits[1].nrows;

    free(dataSplits);

    (rootCopy->sizeLR)[0] = nLeft;
    (rootCopy->sizeLR)[1] = nRight;

    if (nLeft != 0)
    {
        TreeOOBInternal(
            rootCopy->leftChild,
            designMatrixYLeft,
            auxiliaryFeaturesLeft,
            nLeft,
            ncolsDesign,
            ncolsAuxiliary,
            lenOutput);
    }
    else
    {
        FreeTree(rootCopy->leftChild);
        rootCopy->leftChild = NULL;
    }

    if (nRight != 0)
    {
        TreeOOBInternal(
            rootCopy->rightChild,
            designMatrixYRight,
            auxiliaryFeaturesRight,
            nRight,
            ncolsDesign,
            ncolsAuxiliary,
            lenOutput);
    }
    else
    {
        FreeTree(rootCopy->rightChild);
        rootCopy->rightChild = NULL;
    }

    // free the memory
    free(designMatrixYLeft);
    free(designMatrixYRight);
    free(auxiliaryFeaturesLeft);
    free(auxiliaryFeaturesRight);
}

DecisionTreeNode *TreeOOB(
    DecisionTreeNode *root,
    double **designMatrixY,     // the out-of-bag design matrix
    double **auxiliaryFeatures, // the out-of-bag auxiliary features, pseudo-risk estimation should be done in forest step before feeding into this function, similar with GrowTree()
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    size_t lenOutput)
{
    DecisionTreeNode *rootCopy = TreeCopy(root);

    TreeOOBInternal(
        rootCopy,
        designMatrixY,
        auxiliaryFeatures,
        nrows,
        ncolsDesign,
        ncolsAuxiliary,
        lenOutput);

    return rootCopy;
}

// add up all the terminal node splitStat[0] in the decision tree node
double TreeLikelihoodSum(DecisionTreeNode *root)
{
    if (root == NULL)
    {
        return 0.0;
    }

    if (root->flag == 1)
    {
        if (isnan(root->splitStat[0]))
        {
            return 0.0;
        }
        return root->splitStat[0];
    }

    return TreeLikelihoodSum(root->leftChild) + TreeLikelihoodSum(root->rightChild);
}

void FreeTree(DecisionTreeNode *tree)
{
    if (tree == NULL)
    {
        return;
    }
    free(tree->splitStat);
    free(tree->sizeLR);
    free(tree->output);
    FreeTree(tree->leftChild);
    FreeTree(tree->rightChild);
    free(tree);
}

double *TreePredict(
    DecisionTreeNode *root,
    double *CPIU)
{
    if (root->flag == 1)
    {
        return root->output;
    }
    else
    {
        if (CPIU[root->splitIndex] <= root->splitValue)
        {
            return TreePredict(root->leftChild, CPIU);
        }
        else
        {
            return TreePredict(root->rightChild, CPIU);
        }
    }
}

void PrintTreeDot(
    FILE *file,
    DecisionTreeNode *root,
    unsigned int level)
{
    if (root == NULL)
    {
        return;
    }

    if (root->flag == 0)
    {
        fprintf(file, "\tnode%ld [label = \"X%ld <= %.2f \n L(%ld)R(%ld); %.3f \n S(P%.2f;L%.2f;R%.2f)\"] \n", root->nodeId, root->splitIndex, root->splitValue, root->sizeLR[0], root->sizeLR[1], root->splitStat[0], root->splitStat[1], root->splitStat[2], root->splitStat[3]);
        fprintf(file, "\tnode%ld -> node%ld\n", root->nodeId, root->leftChild->nodeId);
        fprintf(file, "\tnode%ld -> node%ld\n", root->nodeId, root->rightChild->nodeId);
    }
    else
    {
        if (root->lenOutput > 1)
        {
            fprintf(file, "\tnode%ld [shape = record label=\"%.2f ...\"]\n", root->nodeId, root->output[0]);
        }
        else
        {
            fprintf(file, "\tnode%ld [shape = record label=\"%.2f\"]\n", root->nodeId, root->output[0]);
        }
    }

    PrintTreeDot(file, root->leftChild, level + 1);
    PrintTreeDot(file, root->rightChild, level + 1);
}

void WriteTreeDotFile(
    DecisionTreeNode *root,
    char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        PRINT_LOCATION();
        printf("Error opening file!");
        exit(1);
    }

    fprintf(file, "digraph {\n");
    fprintf(file, "//Tree ID: %ld\n", root->treeId);
    PrintTreeDot(file, root, 0);
    fprintf(file, "}\n");

    fclose(file);
}

void VimpTree(
    DecisionTreeNode *root,
    double **vimpStat,
    double **vimpFreq)
{
    if (root == NULL)
    {
        return;
    }

    if (root->flag == 0)
    {
        (*vimpStat)[root->splitIndex] += root->splitStat[1];
        (*vimpFreq)[root->splitIndex] += 1;
    }

    VimpTree(root->leftChild, vimpStat, vimpFreq);
    VimpTree(root->rightChild, vimpStat, vimpFreq);
}

void PrintTree(DecisionTreeNode *root)
{
    if (root == NULL)
    {
        return;
    }

    if (root->flag == 0)
    {
        printf("Node %ld, flag is %d, split index %ld, split value %.3f, split stat %.3f, sizeLR %ld %ld\n", root->nodeId, root->flag, root->splitIndex, root->splitValue, root->splitStat[0], root->sizeLR[0], root->sizeLR[1]);
    }
    else
    {
        printf("Node %ld, flag is %d, output %.3f\n", root->nodeId, root->flag, root->output[0]);
    }

    PrintTree(root->leftChild);
    PrintTree(root->rightChild);
}

void SaveTree(
    DecisionTreeNode *root,
    FILE *file)
{
    if (root == NULL)
    {
        return;
    }
    else
    {
        if (root->flag)
        {
            fprintf(file, "%d %ld %ld %ld ", root->flag, root->treeId, root->nodeId, root->lenOutput);
            for (size_t i = 0; i < root->lenOutput; i++)
            {
                fprintf(file, "%.3f ", root->output[i]);
            }
            fprintf(file, "\n");
        }
        else
        {
            fprintf(file, "%d %ld %ld %ld %.3f %.3f %.3f %.3f %.3f %ld %ld\n", root->flag, root->treeId, root->nodeId, root->splitIndex, root->splitValue, root->splitStat[0], root->splitStat[1], root->splitStat[2], root->splitStat[3], root->sizeLR[0], root->sizeLR[1]);
        }
        SaveTree(root->leftChild, file);
        SaveTree(root->rightChild, file);
    }
}

DecisionTreeNode *LoadTree(FILE *file)
{
    int flag;
    long treeId;
    long nodeId;

    if (file == NULL)
    {
        PRINT_LOCATION();
        printf("Error opening file.\n");
        exit(1);
    }

    if (fscanf(file, "%d %ld %ld", &flag, &treeId, &nodeId) != 3)
    {
        fprintf(stderr, "Error reading flag from file.\n");
        return NULL;
    }

    DecisionTreeNode *root = malloc(sizeof(DecisionTreeNode));
    if (root == NULL)
    {
        fprintf(stderr, "Error allocating memory for tree node.\n");
        return NULL;
    }
    root->flag = flag;
    root->treeId = treeId;
    root->nodeId = nodeId;
    root->leftChild = NULL;
    root->rightChild = NULL;
    root->sizeLR = NULL;
    root->output = NULL;
    root->splitStat = NULL;

    if (root->flag)
    {
        if (fscanf(file, "%zu ", &root->lenOutput) != 1)
        {
            PRINT_LOCATION();
            printf("Error: reading leaf node output length from file.\n");
            free(root);
            exit(1);
        }
        root->output = (double *)calloc(root->lenOutput, sizeof(double));
        if (root->output == NULL)
        {
            PRINT_LOCATION();
            printf("Error: allocating memory for leaf node output.\n");
            free(root);
            return NULL;
        }

        for (size_t i = 0; i < (root->lenOutput - 1); i++)
        {
            if (fscanf(file, "%lf ", &root->output[i]) != 1)
            {
                PRINT_LOCATION();
                printf("Error: reading leaf node output from file.\n");
                free(root->output);
                free(root);
                return NULL;
            }
        }
        if (fscanf(file, "%lf \n", &root->output[(root->lenOutput) - 1]) != 1)
        {
            PRINT_LOCATION();
            printf("Error: reading leaf node output from file.\n");
            free(root->output);
            free(root);
            return NULL;
        }
    }
    else
    {
        root->sizeLR = malloc(2 * sizeof(size_t));
        root->splitStat = malloc(4 * sizeof(double));
        if (root->sizeLR == NULL || root->splitStat == NULL)
        {
            PRINT_LOCATION();
            printf("Error: allocating memory for internal node data.\n");
            free(root->sizeLR);
            free(root->splitStat);
            free(root);
            return NULL;
        }
        if (fscanf(file, "%ld %lf %lf %lf %lf %lf %zu %zu\n", &root->splitIndex, &root->splitValue, &root->splitStat[0], &root->splitStat[1], &root->splitStat[2], &root->splitStat[3], &root->sizeLR[0], &root->sizeLR[1]) != 8)
        {
            PRINT_LOCATION();
            printf("Error: reading internal node data from file.\n");
            free(root->sizeLR);
            free(root->splitStat);
            free(root);
            return NULL;
        }
        root->leftChild = LoadTree(file);
        if (root->leftChild == NULL)
        {
            PRINT_LOCATION();
            printf("Error: creating left child node.\n");
            free(root->sizeLR);
            free(root->splitStat);
            free(root);
            return NULL;
        }
        root->rightChild = LoadTree(file);
        if (root->rightChild == NULL)
        {
            PRINT_LOCATION();
            printf("Error: creating right child node.\n");
            free(root->sizeLR);
            free(root->splitStat);
            free(root->leftChild);
            free(root);
            return NULL;
        }
    }

    return root;
}

void PrintAllNodesElementsRecur(
    DecisionTreeNode *root,
    FILE *file,
    NodeElements **nodeElements,
    size_t pathDepth,
    size_t *nthPath)
{
    if (root == NULL)
    {
        return;
    }

    nodeElements[pathDepth]->treeId = root->treeId;
    nodeElements[pathDepth]->nodeId = root->nodeId;
    nodeElements[pathDepth]->flag = root->flag;
    nodeElements[pathDepth]->splitIndex = root->splitIndex;
    nodeElements[pathDepth]->splitValue = root->splitValue;
    nodeElements[pathDepth]->splitStat = root->splitStat;
    nodeElements[pathDepth]->sizeLR = root->sizeLR;
    nodeElements[pathDepth]->output = root->output;
    nodeElements[pathDepth]->lenOutput = root->lenOutput;

    pathDepth++;

    if (root->flag == 1)
    {
        for (int i = 0; i < pathDepth - 1; i++)
        {
            fprintf(
                file,
                "%ld,%ld,%ld,%u,%ld,%lf,%lf,%lf,%lf,%lf,%zu,%zu,\n",
                nodeElements[i]->treeId,
                nodeElements[i]->nodeId,
                *nthPath,
                nodeElements[i]->flag,
                nodeElements[i]->splitIndex,
                nodeElements[i]->splitValue,
                nodeElements[i]->splitStat[0],
                nodeElements[i]->splitStat[1],
                nodeElements[i]->splitStat[2],
                nodeElements[i]->splitStat[3],
                nodeElements[i]->sizeLR[0],
                nodeElements[i]->sizeLR[1]);
        }

        fprintf(
            file,
            "%ld,%ld,%ld,%u,,,,,,,,,",
            nodeElements[pathDepth - 1]->treeId,
            nodeElements[pathDepth - 1]->nodeId,
            *nthPath,
            nodeElements[pathDepth - 1]->flag);

        for (int j = 0; j < nodeElements[pathDepth - 1]->lenOutput - 1; j++)
        {
            fprintf(file, "%lf,", (nodeElements[pathDepth - 1]->output)[j]);
        }

        fprintf(file, "%lf\n", nodeElements[pathDepth - 1]->output[nodeElements[pathDepth - 1]->lenOutput - 1]);

        (*nthPath)++;
    }

    PrintAllNodesElementsRecur(root->leftChild, file, nodeElements, pathDepth, nthPath);
    PrintAllNodesElementsRecur(root->rightChild, file, nodeElements, pathDepth, nthPath);
}

void PrintTreeNodes(
    DecisionTreeNode *root,
    FILE *file)
{
    // Allocate memory for path
    // could print out the tree path up to length of 1000
    NodeElements **nodeElements = malloc(1000 * sizeof(NodeElements *));

    for (int i = 0; i < 1000; i++)
    {
        nodeElements[i] = malloc(sizeof(NodeElements));
    }

    if (nodeElements == NULL)
    {
        PRINT_LOCATION();
        printf("Error: allocating memory for path.\n");
        exit(1);
    }

    size_t nthPath = 0;

    PrintAllNodesElementsRecur(root, file, nodeElements, 0, &nthPath);

    for (int i = 0; i < 1000; i++)
    {
        free(nodeElements[i]);
    }
    free(nodeElements);
}