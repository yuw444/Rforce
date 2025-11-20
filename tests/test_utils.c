#include "unity/unity.h"
#include <stdio.h>
#include <stdlib.h>
#include "../src/utils.h"
#include "../src/pvalue.h"

int *seqInt;
double *seqDouble;
size_t len = 10;
double **array2D;

void SuiteSetup(void)
{
    size_t ncolD = 36;
    seqInt = GetSeqInt(0, ncolD - 1, 1);
    PrintArrayInt(seqInt, 36);
    seqDouble = GetSeqDouble(0, 1.0f, 1.0 / len);

    array2D = Allocate2DArray(2, len);
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < len; j++)
        {
            array2D[i][j] = i * len + j;
        }
    }
}

void SuiteTeardown(void)
{
    free(seqInt);
    free(seqDouble);
    Free2DArray(array2D, 2);
}

/* Unity calls these before/after each test; we don't need per-test setup */
void setUp(void) {}
void tearDown(void) {}

void test_GetSeqInt(void)
{
    TEST_ASSERT_EQUAL_INT(0, seqInt[0]);
    TEST_ASSERT_EQUAL_INT(1, seqInt[1]);
    TEST_ASSERT_EQUAL_INT(2, seqInt[2]);
    TEST_ASSERT_EQUAL_INT(3, seqInt[3]);
    TEST_ASSERT_EQUAL_INT(4, seqInt[4]);
    TEST_ASSERT_EQUAL_INT(5, seqInt[5]);
    TEST_ASSERT_EQUAL_INT(6, seqInt[6]);
    TEST_ASSERT_EQUAL_INT(7, seqInt[7]);
    TEST_ASSERT_EQUAL_INT(8, seqInt[8]);
}

void test_GetSeqDouble(void)
{
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 0.0f, (float)seqDouble[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 0.8f, (float)seqDouble[8]);
}

void test_SampleInt(void)
{
    int *sample1 = SampleInt(seqInt, len, len, 0, 1);
    printf("Sample Int without replacement:\n");
    PrintArrayInt(sample1, len);

    int *temp = SampleInt(seqInt, len, len + 1, 1, 1);
    free(temp);

    int *sample2 = SampleInt(seqInt, len, len, 1, 1);
    printf("Sample Int with replacement:\n");
    PrintArrayInt(sample2, len);

    free(sample1);
    free(sample2);
}

void test_SampleInt2(void)
{
    int *sample2 = GetSeqInt(0, 9, 1);
    int **matrixSample = (int **)malloc(1000 * sizeof(int *));

    for (size_t i = 0; i < 1000; i++)
    {
        matrixSample[i] = SampleInt(sample2, 10, 3, 0, 926);
    }

    FILE *fp = fopen("/home/yu89975/randomForest/test/matrixSample.csv", "w");
    if (fp)
    {
        WriteCSVInt(matrixSample, fp, 1000, 3);
        fclose(fp);
    }

    for (size_t i = 0; i < 1000; i++) free(matrixSample[i]);
    free(matrixSample);
    free(sample2);
}

void test_SampleDouble(void)
{
    double *sample3 = SampleDouble(seqDouble, len, len, 0, 1);
    printf("Sample Double without replacement:\n");
    PrintArrayDouble(sample3, len);

    double *sample4 = SampleDouble(seqDouble, len, len, 1, 1);
    printf("Sample Double with replacement:\n");
    PrintArrayDouble(sample4, len);

    double *sample5 = SampleDouble(seqDouble, len, len + 1, 0, 1);
    TEST_ASSERT_NULL(sample5);

    free(sample3);
    free(sample4);
}

void test_Unique(void)
{
    printf("Unique test:\n");
    PrintArrayDouble(seqDouble, len);
    ElementStruct *unique = Unique(seqDouble, len);
    TEST_ASSERT_EQUAL_INT((int)len, unique->nElements);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 0.0f, (float)unique->Elements[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 0.2f, (float)unique->Elements[2]);
    free(unique->Elements);
    free(unique);

    double *tt1 = GetSeqDouble(1.0, 10.0, 1.0);
    tt1[2] = NA_DOUBLE;
    tt1[7] = NA_DOUBLE;

    printf("Unique is handling NA_DOUBLE:\n\nBefore:\n");
    ElementStruct *unique2 = Unique(tt1, 10);
    PrintArrayDouble(unique2->Elements, unique2->nElements);

    free(tt1);
    free(unique2->Elements);
    free(unique2);
}

void test_Quantile(void)
{
    double quantils[] = {0.0, 0.25, 0.5, 0.75, 1.0};
    double *quantile = Quantile(seqDouble, len, quantils, 5);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 0.0f, (float)quantile[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 0.2f, (float)quantile[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 0.4f, (float)quantile[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 0.6f, (float)quantile[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 0.9f, (float)quantile[4]);
    free(quantile);
}

void test_SplitCandidates(void)
{
    PrintArrayDouble(seqDouble, len);
    ElementStruct *splitCandidates = SplitCandidates(seqDouble, len, 3);
    PrintArrayDouble(splitCandidates->Elements, splitCandidates->nElements);
    TEST_ASSERT_EQUAL_INT(3, splitCandidates->nElements);
    free(splitCandidates->Elements);
    free(splitCandidates);
}

void test_GetCol(void)
{
    double *col = GetCol(array2D, 2, len, 1);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 1.0f, (float)col[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 11.0f, (float)col[1]);
    free(col);
}

void test_GetNrowCSV(void)
{
    char *filename = "/home/yu89975/randomForest/test/test.csv";
    int nrow = GetNrowCSV(filename);
    TEST_ASSERT_EQUAL_INT(4, nrow);
}

void test_GetNcolCSV(void)
{
    char *filename = "/home/yu89975/randomForest/test/test.csv";
    int ncol = GetNcolCSV(filename);
    TEST_ASSERT_EQUAL_INT(3, ncol);
}

void test_ReadCSV(void)
{
    char *filename = "/home/yu89975/randomForest/test/test.csv";
    double **data = Allocate2DArray(4, 3);
    ReadCSV(filename, &data, 0);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 1.0f, (float)data[0][0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 2.0f, (float)data[0][1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 3.0f, (float)data[0][2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 4.0f, (float)data[1][0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 5.0f, (float)data[1][1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 6.0f, (float)data[1][2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 7.0f, (float)data[2][0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 8.0f, (float)data[2][1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 9.0f, (float)data[2][2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 10.0f, (float)data[3][0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 11.0f, (float)data[3][1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 12.0f, (float)data[3][2]);
    Free2DArray(data, 4);

    char *filename1 = "/home/yu89975/randomForest/data/testNA.csv";
    double **dataNA = Allocate2DArray(2, 12);
    ReadCSV(filename1, &dataNA, 0);
    PrintArrayDouble(dataNA[0], 12);
    TEST_ASSERT_TRUE(isnan(dataNA[0][2]));
    Free2DArray(dataNA, 2);
}

void test_WriteCSV(void)
{
    char *filename = "/home/yu89975/randomForest/test/test.csv";
    double **data = Allocate2DArray(4, 3);
    ReadCSV(filename, &data, 0);
    FILE *fp = fopen("/home/yu89975/randomForest/test/test2.csv", "w");
    if (fp)
    {
        WriteCSV(data, fp, 4, 3);
        fclose(fp);
    }
    double **data2 = Allocate2DArray(4, 3);
    ReadCSV("/home/yu89975/randomForest/test/test2.csv", &data2, 0);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 1.0f, (float)data2[0][0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 2.0f, (float)data2[0][1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 3.0f, (float)data2[0][2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 4.0f, (float)data2[1][0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 5.0f, (float)data2[1][1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 6.0f, (float)data2[1][2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 7.0f, (float)data2[2][0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 8.0f, (float)data2[2][1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 9.0f, (float)data2[2][2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 10.0f, (float)data2[3][0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 11.0f, (float)data2[3][1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 12.0f, (float)data2[3][2]);
    Free2DArray(data, 4);
    Free2DArray(data2, 4);
}

void test_Factorial(void)
{
    int n = 22;
    double fact = lgamma(n + 1);
    printf("test are %.3f\n", fact);
}

void test_ThreeDArray(void)
{
    double ***array = Allocate3DArray(2, 3, 1);
    Free3DArray(array, 2, 3);
}

void test_NthColMean(void)
{
    double **array = Allocate2DArray(10, 10);
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            array[i][j] = 1.0;
        }
    }

    printf("Array for NthColMean test:\n");

    array[0][1] = NA_DOUBLE;
    array[0][9] = NA_DOUBLE;
    array[1][1] = 10;
    array[1][9] = 10;

    for (int i = 0; i < 10; i++)
    {
        PrintArrayDouble(array[i], 10);
    }

    double meanFirstCol = NthColMean(array, 10, 0);
    double meanLastCol = NthColMean(array, 10, 9);
    double meanSecondCol = NthColMean(array, 10, 1);

    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 1.0f, (float)meanFirstCol);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 2.0f, (float)meanLastCol);
    TEST_ASSERT_FLOAT_WITHIN(1e-6F, 2.0f, (float)meanSecondCol);

    Free2DArray(array, 10);
}

void test_ReadingNA(void)
{
    char *designMatrixYFile = "/home/yu89975/randomForest/data/testNA.csv";
    int nrowD = GetNrowCSV(designMatrixYFile);
    int ncolD = GetNcolCSV(designMatrixYFile);
    double **designMatrixY = Allocate2DArray(nrowD, ncolD);
    ReadCSV(designMatrixYFile, &designMatrixY, 0);

    int *nNA = (int *)calloc(ncolD, sizeof(int));
    for (int i = 0; i < nrowD; i++)
    {
        for (int j = 0; j < ncolD; j++)
        {
            if (isnan(designMatrixY[i][j]))
            {
                nNA[j]++;
            }
        }
    }

    PrintArrayInt(nNA, ncolD);

    free(nNA);
    Free2DArray(designMatrixY, nrowD);
}

void test_pvalueAdjust(void)
{
    double pvalues[] = {0.01, 0.08, 0.03, 0.04, 0.05};
    double *pvaluesAdj = PAdjust(pvalues, 5, "BH");
    PrintArrayDouble(pvaluesAdj, 5);
    free(pvaluesAdj);
}

void test_Permutate(void)
{
    double x[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double *xperm = Permute(x, 5, 1);
    PrintArrayDouble(xperm, 5);
    free(xperm);
}

void test_ColsPermute(void)
{
    double **x = Allocate2DArray(10, 10);
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            x[i][j] = i * 3 + j + 1;
        }
    }

    unsigned int *colsToPermute = GetSeqInt(3, 5, 1);

    double **xperm = ColsPermute(x, 10, 10, colsToPermute, 3, 926);

    for (int i = 0; i < 10; i++)
    {
        PrintArrayDouble(x[i], 10);
    }

    for (int i = 0; i < 10; i++)
    {
        PrintArrayDouble(xperm[i], 10);
    }

    Free2DArray(x, 10);
    Free2DArray(xperm, 10);
    free(colsToPermute);
}


int main(void)
{
    SuiteSetup();
    UNITY_BEGIN();

    RUN_TEST(test_GetSeqInt);
    RUN_TEST(test_GetSeqDouble);
    RUN_TEST(test_SampleInt);
    RUN_TEST(test_SampleInt2);
    RUN_TEST(test_SampleDouble);
    RUN_TEST(test_Unique);
    RUN_TEST(test_Quantile);
    RUN_TEST(test_SplitCandidates);
    RUN_TEST(test_GetCol);
    RUN_TEST(test_GetNrowCSV);
    RUN_TEST(test_GetNcolCSV);
    RUN_TEST(test_ReadCSV);
    RUN_TEST(test_WriteCSV);
    RUN_TEST(test_Factorial);
    RUN_TEST(test_ThreeDArray);
    RUN_TEST(test_NthColMean);
    RUN_TEST(test_ReadingNA);
    RUN_TEST(test_pvalueAdjust);
    RUN_TEST(test_Permutate);
    RUN_TEST(test_ColsPermute);

    int result = UNITY_END();
    SuiteTeardown();
    return result;
}