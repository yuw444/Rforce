#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include "argparse.h"
#include "forest.h"

#define ARRAY_SIZE(x) (sizeof(x) / sizeof(x[0]))

static const char *const usages[] = {
    "Rforce [subcommands] <options> \nsubcommands: \n\ttrain: train a composite endpoint forest \n\tpredict: predict according to a composite endpoint forest and observations",
    NULL,
};

struct cmd_struct
{
    const char *cmd;
    int (*fn)(int, const char **);
};

int cmd_train(int argc, const char **argv)
{
    char *path_designmatrixY = NULL;
    char *path_auxiliary = NULL;
    char *path_unitsOfCPIU = NULL;
    char *path_out = ".";
    char *pathVarIds = NULL;
    _padjust = "BH";

    size_t lenOutPut;
    size_t maxDepth = 10;
    size_t minNodeSize = 5;
    float minGain = 0.f;
    size_t mtry = 0;
    size_t nsplits = 10;
    size_t nTrees = 200;
    unsigned int seed = 926;
    size_t n_threads = 8;
    size_t nVars = 0;
    unsigned int *varIds = NULL;

    struct argparse_option options[] = {
        OPT_HELP(),
        // input and ouput
        OPT_STRING('d', "designMatrixY", &path_designmatrixY, "required; path to design matrix", NULL, 0, 0),
        OPT_STRING('a', "auxiliary", &path_auxiliary, "required; path to auxiliary features", NULL, 0, 0),
        OPT_STRING('u', "unitsOfCPIU", &path_unitsOfCPIU, "required; path to unitsOfCPIU", NULL, 0, 0),
        OPT_STRING('o', "out", &path_out, "path to output, optional, default is current working directory", NULL, 0, 0),
        OPT_INTEGER('v', "verbose", &_verbose, "level of verbose, 0, 1, 2, 3, default 0", NULL, 0, 0),
        // model parameters
        OPT_INTEGER('m', "maxDepth", &maxDepth, "max depth of tree; optional, default is 10", NULL, 0, 0),
        OPT_INTEGER('n', "minNodeSize", &minNodeSize, "min node size; optional, default is 2 * len(unitsOfCPIU) - 1", NULL, 0, 0),
        OPT_FLOAT('g', "gain", &minGain, "min gain; default is 0.0 when using likelihood based split, -log10(p-value = 0.05) = 1.3 when using GEE and Aysmpotic", NULL, 0, 0),
        OPT_INTEGER('t', "mtry", &mtry, "number of variables to try during splitting; optional, default is the sqrt of number of variables", NULL, 0, 0),
        OPT_INTEGER('s', "nsplits", &nsplits, "number of splitting to try at each variable; optional, default is 10", NULL, 0, 0),
        OPT_INTEGER('r', "nTrees", &nTrees, "nTrees; optional, default is 200", NULL, 0, 0),
        OPT_INTEGER('e', "seed", &seed, "seed for random number generator; optional, deafult is 926", NULL, 0, 0),
        OPT_INTEGER('p', "nPerms", &_nPerms, "number of permutations for permutation importance for each variable and each tree; optional, default is 10", NULL, 0, 0),

        // optional
        OPT_INTEGER('u', "nVars", &nVars, "number of variables in the design matrix; optional, default is the number of columns of variables", NULL, 0, 0),
        OPT_STRING('i', "pathVarIds", &pathVarIds, "a vector of variable IDs, from 0 to (nVars -1). the consecutive repeated Ids indicates a single categorial variable when calculating vimpPermuted", NULL, 0, 0),
        
        // node prediction
        OPT_INTEGER('k', "k", &_k, "parameter for bayesian estimator in leaf node output, default is 4, bigger means less info borrow", NULL, 0, 0),

        /**********************Split Rules************************/
        // Original RF-SLAM
        OPT_BOOLEAN('L', "long", &_longformat, "one time to event in one row, one patient may take multiple rows, which corresponds to the original RF-SLAM", NULL, 0, 0),

        // Quasi-likelihood approaches
        OPT_BOOLEAN('N', "nopseudo", &_noPseudo, "do not do pseudo risk time estimation, keep the original observed risk time in the auxiliary matrix", NULL, 0, 0),
        OPT_BOOLEAN('P', "pseudorisk1", &_pseudoRisk1, "use the original pseudo-risk time in the auxiliary matrix that calculated at population level", NULL, 0, 0),
        OPT_BOOLEAN('B', "pseudorisk2", &_pseudoRisk2, "default, re-calculate pseudo risk time at each tree level", NULL, 0, 0),
        OPT_BOOLEAN('D', "dynamicrisk", &_dynamicRisk, "use the dynamic pseudo-risk time estimation at each split", NULL, 0, 0),
        
        OPT_BOOLEAN('F', "nophi", &_noPhi, "do not do the phi estimation, use the fixed phi = 1", NULL, 0, 0),
        OPT_BOOLEAN('P', "phi1", &_phi1, "phi is calculated at the population level", NULL, 0, 0),
        OPT_BOOLEAN('H', "phi2", &_phi2, "default, calculate phi at each tree level", NULL, 0, 0),
        OPT_BOOLEAN('Y', "dynamicphi", &_dynamicPhi, "use the dynamic phi estimation at each split", NULL, 0, 0),

        // GEE approaches
        OPT_BOOLEAN('G', "gee", &_gee, "use the GEE approach", NULL, 0, 0),
        OPT_STRING('A', "padjust", &_padjust, "p-value adjustment method,'bonferroni', 'holm', 'hochberg', 'hommel', 'BH', 'BY', 'none'; default is BH,", NULL, 0, 0),
        OPT_BOOLEAN('I', "interaction", &_interaction, "add interaction term for GEE, default is NULL", NULL, 0, 0),
        // Asympotic approaches
        OPT_BOOLEAN('S', "asym", &_asympotic, "use the asympotic approach", NULL, 0, 0),

        // OpenMP
        OPT_INTEGER('T', "threads", &n_threads, "number of threads for parallel computing; optional, default is 8", NULL, 0, 0),
        OPT_END()
    };

    if (pathVarIds != NULL)
    {
        nVars = GetNrowCSV(pathVarIds);
        varIds = (unsigned int *)calloc(nVars, sizeof(unsigned int));
        double **varIds_temp = Allocate2DArray(nVars, 1);
        ReadCSV(pathVarIds, &varIds_temp, 0);
        for (size_t i = 0; i < nVars; i++)
        {
            varIds[i] = (unsigned int)varIds_temp[i][0];
        }
        Free2DArray(varIds_temp, nVars);
    }
    
    struct argparse argparse;
    argparse_init(&argparse, options, usages, 0);

    argc = argparse_parse(&argparse, argc, argv);

    //set up the number of threads
    omp_set_num_threads(n_threads);

    // check if inputs are compatiable
    int sum_risk = _noPseudo + _pseudoRisk1 + _pseudoRisk2 + _dynamicRisk;
    int sum_phi = _noPhi + _phi1 + _phi2 + _dynamicPhi;

    if (_longformat == 0 && _gee == 0 && (sum_risk != 1 || sum_phi != 1))
    {
        printf("Error: only one of risk options --nopseudo, --pseudorisk1, --pseudorisk2, --dynamicrisk can be used.\n");
        printf("Error: only one of phi options --nophi, --phi1, --phi2, --dynamicphi can be used.\n");
        exit(1);
    }

    if (_gee == 1 && (sum_risk != 1 || sum_phi != 1))
    {
        printf("Error: only one of risk options --nopseudo, --pseudorisk1, --pseudorisk2, --dynamicrisk can be used.\n");
        printf("Error: no phi-related can be used here.\n");
        exit(1);
    }

    if (_longformat == 1 && (_gee == 1 || _asympotic == 1))
    {
        printf("Warning: --longformat is not compatiable with other arguments, only --longformat is used.\n");
    }

    const char p_adjust_method[7][20] = {"bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "none"};
    
    // check if _padjust is valid
    int valid_p_adjust = 0;
    for (int i = 0; i < 7; i++)
    {
        if (strcmp(_padjust, p_adjust_method[i]) == 0)
        {
            valid_p_adjust = 1;
            break;
        }
    }

    if (valid_p_adjust == 0)
    {
        printf("Error: invalid p-value adjustment method, please use one of 'bonferroni', 'holm', 'hochberg', 'hommel', 'BH', 'BY', 'none'\n");
        exit(1);
    }

    // check if all required arguments are provided
    if (path_designmatrixY == NULL || path_auxiliary == NULL || path_out == NULL)
    {
        printf("Error: required arguments are missing\n");
        argparse_usage(&argparse);
        return 0;
    }

    // check if path_out exists
    if (access(path_out, F_OK) == -1)
    {
        printf("Error: path_out does not exist\n");
        exit(1);
    }
    else
    {
        char path_out_temp[1024];
        sprintf(
            path_out_temp,
            "%s/%s_tree%zu_depth%zu_mtry%zu_Gain%.3f_nsplits%zu_nodesize%zu_seed%d",
            path_out,
            "forest",
            nTrees,
            maxDepth,
            mtry,
            minGain,
            nsplits,
            minNodeSize,
            seed);

        if (access(path_out_temp, F_OK) == 0)
        {
            printf("Warning: %s already exists, please move it other place.\n", path_out_temp);
            printf("Options:\n");
            printf("0. Exit.\n");
            printf("1. Use another path for output.\n");
            printf("2. Delete the existing folder %s.\n", path_out_temp);

            int choice;
            printf("Please enter your choice: ");
            scanf("%d", &choice);

            if (choice == 1)
            {
                printf("Please enter the path for output: ");
                scanf("%s", path_out);
            }
            else if (choice == 2)
            {
                printf("Deleting %s...\n", path_out_temp);
                if (DeleteDir(path_out_temp) == -1)
                {
                    printf("Error: failed to delete %s\n", path_out_temp);
                    exit(1);
                }
            }
            else if (choice == 0)
            {
                printf("Exiting...\n");
                exit(0);
            }
            else
            {
                printf("Invalid choice.\n");
                exit(1);
            }
        }
    }

    // read designMatrixY
    size_t nrowD = GetNrowCSV(path_designmatrixY);
    size_t ncolD = GetNcolCSV(path_designmatrixY);
    mtry = mtry == 0 ? (size_t)sqrt(ncolD) : mtry;
    double **designMatrixY = Allocate2DArray(nrowD, ncolD);
    ReadCSV(path_designmatrixY, &designMatrixY, 0);

    // read auxiliary
    size_t nrowA = GetNrowCSV(path_auxiliary);
    size_t ncolA = GetNcolCSV(path_auxiliary);
    double **auxiliary = Allocate2DArray(nrowA, ncolA);
    ReadCSV(path_auxiliary, &auxiliary, 0);

    // read unitsOfCPIU
    size_t nUnits;
    double *unitsOfCPIU;
    size_t nrowU = GetNrowCSV(path_unitsOfCPIU);
    size_t ncolU = GetNcolCSV(path_unitsOfCPIU);
    double **unitsOfCPIU_temp = Allocate2DArray(nrowU, ncolU);

    ReadCSV(path_unitsOfCPIU, &unitsOfCPIU_temp, 0);

    if (ncolU != 1 && nrowU != 1)
    {
        printf("%s should have only one column or one row\n", path_unitsOfCPIU);
        Free2DArray(unitsOfCPIU_temp, nrowU);
        exit(1);
    }

    // unitsOfCPIU = GetCol(unitsOfCPIU_temp, nrowU, ncolU, ncolU-1);

    if (ncolU == 1)
    {
        unitsOfCPIU = GetCol(unitsOfCPIU_temp, nrowU, ncolU, ncolU - 1);
        PrintArrayDouble(unitsOfCPIU, nrowU);
        nUnits = nrowU;
    }
    else
    {
        unitsOfCPIU = (double *)calloc(ncolU, sizeof(double));
        memcpy(unitsOfCPIU, unitsOfCPIU_temp[0], ncolU * sizeof(double));
        PrintArrayDouble(unitsOfCPIU, ncolU);
        nUnits = ncolU;
    }

    MATRIX *_Rin = VC_GEE_create_matrix(nUnits, nUnits, PERMANENT);

    Free2DArray(unitsOfCPIU_temp, nrowU);

    // guard for gee model, may be unnecessary
    // if (minNodeSize == 0 || (_gee == 1 && minNodeSize < 2 * nUnits - 1))
    // {
    //     minNodeSize = 2 * nUnits - 1;
    // }

    if (_verbose >= 1)
    {
        printf("\n\n##############################\n\n");
        printf("Path to design matrix: %s\n", path_designmatrixY);
        printf("Path to auxiliary: %s\n", path_auxiliary);
        printf("Path to unitsOfCPIU: %s\n", path_unitsOfCPIU);
        printf("Path to output: %s\n", path_out);

        printf("\n\n##############################\n\n");
        printf("The dimensions of design matrix is : %zu %zu\n", nrowD, ncolD);
        printf("The dimensions of auxiliary is : %zu %zu\n", nrowA, ncolA);
        printf("The number of CPIU units is : %zu\n", nUnits);
        printf("The units of CPIU are : ");
        for (size_t i = 0; i < nUnits; i++)
        {
            printf("%f ", unitsOfCPIU[i]);
        }
        printf("\n");

        printf("\n\n##############################\n\n");
        printf("The length of output at leaf node is : %zu\n", nUnits);
        printf("The max depth of tree is : %zu\n", maxDepth);
        printf("The min node size is : %zu\n", minNodeSize);
        printf("The min gain is : %f\n", minGain);
        printf("The number of variables to try during splitting is : %zu\n", mtry);
        printf("The number of splitting to try at each variable is : %zu\n", nsplits);
        printf("The number of trees is : %zu\n", nTrees);
        printf("The seeds for random composite endpoint forest is: %ld\n", seed);
    }
    // random composite endpoint forest
    printf("\n\n##############################\n\n");
    printf("Training the random composite endpoint forest, please wait...\n");

    RandomSurvivalForest *rsf = NULL;
    _treePhi = Allocate2DArray(nTrees, nUnits);

    if (_longformat == 0)
    {
        ncolD = ncolD - nUnits;
        lenOutPut = nUnits;

        // if nVars is not provided, set it to the number of columns of design matrix
        if (nVars == 0)
        {
            nVars = ncolD;
            varIds = GetSeqInt(0, nVars - 1, 1);
        }

        if (_asympotic || _gee)
        {
            if (_asympotic)
            {
                rsf = RandomForest(
                    AsympoticDiffTestStatsNew,
                    LeafOutputInterval,
                    designMatrixY,
                    auxiliary,
                    nrowD,
                    ncolD,
                    ncolA,
                    varIds,
                    nVars,
                    unitsOfCPIU,
                    nUnits,
                    lenOutPut,
                    maxDepth,
                    minNodeSize,
                    minGain,
                    mtry,
                    nsplits,
                    nTrees,
                    seed);
            }
            else
            {
                rsf = RandomForest(
                    GEERule,
                    LeafOutputInterval,
                    designMatrixY,
                    auxiliary,
                    nrowD,
                    ncolD,
                    ncolA,
                    varIds,
                    nVars,
                    unitsOfCPIU,
                    nUnits,
                    lenOutPut,
                    maxDepth,
                    minNodeSize,
                    minGain,
                    mtry,
                    nsplits,
                    nTrees,
                    seed);
            }
        }
        else
        {
            rsf = RandomForest(
                QuasiPoissonLikelihood,
                LeafOutputInterval,
                designMatrixY,
                auxiliary,
                nrowD,
                ncolD,
                ncolA,
                varIds,
                nVars,
                unitsOfCPIU,
                nUnits,
                lenOutPut,
                maxDepth,
                minNodeSize,
                minGain,
                mtry,
                nsplits,
                nTrees,
                seed);
        }
    }
    else
    {
        // rfslam approach, not need pseudo risk time estimation
        ncolD = ncolD - 1;
        lenOutPut = 1;

        // if nVars is not provided, set it to the number of columns of design matrix
        if (nVars == 0)
        {
            nVars = ncolD;
            varIds = GetSeqInt(0, nVars - 1, 1);
        }

        rsf = RandomForest(
            PoissonLikelihood,
            LeafOutput,
            designMatrixY,
            auxiliary,
            nrowD,
            ncolD,
            ncolA,
            varIds,
            nVars,
            unitsOfCPIU,
            nUnits,
            lenOutPut,
            maxDepth,
            minNodeSize,
            minGain,
            mtry,
            nsplits,
            nTrees,
            seed);
    }

    // save the model
    printf("\n\n##############################\n\n");
    printf("Saving the model to %s, please wait...\n", path_out);
    SaveSurvivalForest(rsf, path_out);

    // print the paths of the forest
    printf("\n\n##############################\n\n");
    printf("Saving the paths of the forest to %s/paths.txt, please wait...\n", path_out);
    char *path_paths = malloc(1024 * sizeof(char));
    sprintf(path_paths, "%s/%s_tree%zu_depth%zu_mtry%zu_Gain%.3f_nsplits%zu_nodesize%zu_seed%d/paths.txt",
            path_out,
            "forest",
            rsf->nTrees,
            rsf->maxDepth,
            rsf->mtry,
            rsf->minGain,
            rsf->nsplits,
            rsf->minNodeSize,
            rsf->seed);
    FILE *fp = fopen(path_paths, "w");
    PrintForestPaths(rsf, fp);
    fclose(fp);
    free(path_paths);

    // free memory
    FreeSurvivalForest(rsf);

    Free2DArray(designMatrixY, nrowD);
    Free2DArray(auxiliary, nrowA);
    Free2DArray(_treePhi, nTrees);
    VC_GEE_destroy_matrix(&_Rin);

    free(unitsOfCPIU);

    printf("Saving completed.\n");

    return 0;
}

int cmd_predict(int argc, const char **argv)
{
    char *path_model = NULL;
    char *path_test = NULL;
    char *path_out = ".";

    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_STRING('m', "model", &path_model, "required; path to model", NULL, 0, 0),
        OPT_STRING('t', "test", &path_test, "required; path to test data", NULL, 0, 0),
        OPT_STRING('o', "out", &path_out, "path to output; optional, default is current working directory", NULL, 0, 0),

        OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usages, 0);
    argc = argparse_parse(&argparse, argc, argv);

    // check if all required arguments are provided
    if (path_model == NULL || path_test == NULL)
    {
        return 1;
    }

    // read model
    printf("\n\n##############################\n\n");
    printf("Loading the model from %s, please wait...\n\n", path_model);

    if (access(path_model, F_OK) == -1)
    {
        printf("Error: %s does not exist\n", path_model);
        exit(1);
    }

    RandomSurvivalForest *rsf = LoadSurvivalForest(path_model);

    // read test data
    size_t ncolT = GetNcolCSV(path_test);
    if (ncolT < rsf->nVars)
    {
        printf("Error: the number of columns in test data is less than the number of variables in training data\n");
        return 1;
    }
    size_t nrowT = GetNrowCSV(path_test);
    double **test = Allocate2DArray(nrowT, ncolT);

    printf("\n\n##############################\n\n");
    printf("Loading the test data from %s, please wait...\n\n", path_test);
    ReadCSV(path_test, &test, 0);

    // predict
    double **predict = SurvivalForestPredict(
        rsf,
        test,
        nrowT);

    // save the prediction
    char *path_predict = (char *)malloc(sizeof(char) * 300);
    sprintf(path_predict, "%s/test_prediction.csv", path_out);
    FILE *fp = fopen(path_predict, "wb");
    WriteCSV(predict, fp, nrowT, rsf->lenOutput);

    printf("Prediction saved to %s\n", path_predict);

    // free memory
    FreeSurvivalForest(rsf);
    Free2DArray(test, nrowT);
    Free2DArray(predict, nrowT);
    free(path_predict);
    fclose(fp);

    return 0;
}

static struct cmd_struct commands[] = {
    {"train", cmd_train},
    {"predict", cmd_predict},
};

int main(int argc, const char **argv)
{
    clock_t startTime = clock();

    struct argparse argparse;
    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_END(),
    };
    argparse_init(&argparse, options, usages, ARGPARSE_STOP_AT_NON_OPTION);
    argc = argparse_parse(&argparse, argc, argv);
    if (argc < 1)
    {
        argparse_usage(&argparse);
        return 1;
    }

    /* Try to run command with args provided. */
    struct cmd_struct *cmd = NULL;
    for (int i = 0; i < ARRAY_SIZE(commands); i++)
    {
        if (!strcmp(commands[i].cmd, argv[0]))
        {
            cmd = &commands[i];
        }
    }
    if (cmd)
    {
        return cmd->fn(argc, argv);
    }

    clock_t endTime = clock();

    printf("Elapsed time: %f seconds\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);

    return 0;
}
