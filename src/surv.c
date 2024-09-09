
#include "surv.h"

void FreeKMResult(KMResult *result)
{
    if (result == NULL)
    {
        return;
    }
    free(result->uniqueTime);
    free(result->survivalProb);
    free(result);
}

// fit a Kaplan-Meier model
KMResult *KaplanMeier(
    double *time,
    double *status,
    int nObs)
{
    // get the time points where the event happened
    // +1 is for padding 0 in the front
    double *timeEvent = (double *)malloc(nObs * sizeof(double));
    if (timeEvent == NULL)
    {
        PRINT_LOCATION();
        printf("Error: failed to allocate memory for timeEvent!\n");
        exit(1);
    }
    int nEvent = 0;

    ElementStruct *uniqueStatus = Unique(status, nObs);
    if (uniqueStatus->nElements > 2)
    {
        PRINT_LOCATION();
        printf("Error: status should only have at most two unique values.\n The current unique values are:\n");
        PrintArrayDouble(uniqueStatus->Elements, uniqueStatus->nElements);
        exit(1);
    }
    free(uniqueStatus->Elements);
    free(uniqueStatus);

    for (int i = 0; i < nObs; i++)
    {
        if (fabs(status[i] - 1) < 1e-6)
        {
            timeEvent[nEvent] = time[i];
            nEvent++;
        }
    }

    if (nEvent == 0)
    {
        double maxTime = Max(time, nObs);
        size_t nTime = 2;
        double *uniqueTime = (double *)malloc(nTime * sizeof(double));
        uniqueTime[0] = 0;
        uniqueTime[1] = maxTime;
        double *survivalProb = (double *)malloc(nTime * sizeof(double));
        survivalProb[0] = 1;
        survivalProb[1] = 1;

        KMResult *result = (KMResult *)malloc(sizeof(KMResult));
        result->nTime = nTime;
        result->uniqueTime = uniqueTime;
        result->survivalProb = survivalProb;
        free(timeEvent);
        return result;
    }

    // printf("nEvent: %d\n", nEvent);

    // get unique time points and their length in sorted order of event time
    ElementStruct *uniqueTime = Unique(timeEvent, nEvent);

    // get the number of at risk and events at each event time point
    int *nAtRisk = (int *)calloc(uniqueTime->nElements, sizeof(int));
    int *nEvents = (int *)calloc(uniqueTime->nElements, sizeof(int));

    for (int i = 0; i < nObs; i++)
    {
        for (int j = 0; j < uniqueTime->nElements; j++)
        {
            if (time[i] >= uniqueTime->Elements[j])
            {
                nAtRisk[j]++;
            }
            if (fabs(time[i] - uniqueTime->Elements[j]) < 1e-8 && status[i] == 1)
            {
                nEvents[j]++;
            }
        }
    }

    // get the Kaplan-Meier estimate
    double *kmEstimate = (double *)malloc(uniqueTime->nElements * sizeof(double));
    if (kmEstimate == NULL)
    {
        PRINT_LOCATION();
        printf("Error: failed to allocate memory for kmEstimate!\n");
        exit(1);
    }
    kmEstimate[0] = 1 - (double)nEvents[0] / nAtRisk[0];

    for (int i = 1; i < uniqueTime->nElements; i++)
    {
        kmEstimate[i] = kmEstimate[i - 1] * (1 - (double)nEvents[i] / nAtRisk[i]);
    }

    // Allocate memory for the result struct
    KMResult *result = (KMResult *)malloc(sizeof(KMResult));
    if (result == NULL)
    {
        PRINT_LOCATION();
        printf("Error: failed to allocate memory for result!\n");
        exit(1);
    }

    result->nTime = uniqueTime->nElements;
    result->uniqueTime = uniqueTime->Elements;
    result->survivalProb = kmEstimate;

    // free memory
    free(timeEvent);
    free(nAtRisk);
    free(nEvents);
    free(uniqueTime);
    return result;
}

double StepFunction(double *x, double *y, size_t len, double x0)
{
    if (x0 < x[0])
    {
        return y[0];
    }

    int i;
    for (i = 0; i < len; i++)
    {
        if (x0 < x[i])
        {
            break;
        }
    }
    return y[i - 1];
}

KMResult *Wt(
    KMResult *Gt,
    double time,
    double status)
{
    if (fabs(status) < 1e-9)
    {
        // printf("status is 0\n");
        KMResult *result = (KMResult *)malloc(sizeof(KMResult));
        if (result == NULL)
        {
            PRINT_LOCATION();
            printf("Error: failed to allocate memory for result!\n");
            exit(1);
        }
        result->nTime = 2;
        result->uniqueTime = (double *)malloc(2 * sizeof(double));
        result->survivalProb = (double *)malloc(2 * sizeof(double));

        result->uniqueTime[0] = 0;
        result->uniqueTime[1] = time;

        result->survivalProb[0] = 1;
        result->survivalProb[1] = 0;

        return result;
    }
    // printf("status is 1 and less than the first time point\n");
    if (time < Gt->uniqueTime[0] || fabs(time - Gt->uniqueTime[0]) < 1e-9)
    {
        // printf("time is less than the first unique time\n");
        KMResult *result = (KMResult *)malloc(sizeof(KMResult));
        if (result == NULL)
        {
            PRINT_LOCATION();
            printf("Error: failed to allocate memory for result!\n");
            exit(1);
        }
        result->nTime = Gt->nTime + 1;
        result->uniqueTime = (double *)malloc(result->nTime * sizeof(double));
        result->survivalProb = (double *)malloc(result->nTime * sizeof(double));

        memcpy(result->uniqueTime + 1, Gt->uniqueTime, Gt->nTime * sizeof(double));
        memcpy(result->survivalProb + 1, Gt->survivalProb, Gt->nTime * sizeof(double));

        result->uniqueTime[0] = 0;
        result->survivalProb[0] = 1;

        return result;
    }
    // printf("status is 1 and greater than the last time point\n");
    if (time > Gt->uniqueTime[Gt->nTime - 1] || fabs(time - Gt->uniqueTime[Gt->nTime - 1]) < 1e-9)
    {
        // printf("time is greater than the last unique time\n");
        KMResult *result = (KMResult *)malloc(sizeof(KMResult));
        if (result == NULL)
        {
            PRINT_LOCATION();
            printf("Error: failed to allocate memory for result!\n");
            exit(1);
        }
        result->nTime = 2;
        result->uniqueTime = (double *)malloc(2 * sizeof(double));
        result->survivalProb = (double *)malloc(2 * sizeof(double));

        result->uniqueTime[0] = 0;
        result->uniqueTime[1] = time;

        result->survivalProb[0] = 1;
        result->survivalProb[1] = 1;

        return result;
    }

    // printf("time is in between\n");
    size_t nGeTime = 0;
    for (int i = 0; i < Gt->nTime; i++)
    {
        if (Gt->uniqueTime[i] < time)
        {
            nGeTime++;
        }
    }

    size_t nNeeds = Gt->nTime - nGeTime + 1;
    KMResult *result = (KMResult *)malloc(sizeof(KMResult));
    if (result == NULL)
    {
        PRINT_LOCATION();
        printf("Error: failed to allocate memory for result!\n");
        exit(1);
    }
    result->nTime = nNeeds;
    result->uniqueTime = (double *)malloc(nNeeds * sizeof(double));
    result->survivalProb = (double *)malloc(nNeeds * sizeof(double));
    result->uniqueTime[0] = 0;
    result->survivalProb[0] = 1;

    for (int i = 1; i < nNeeds; i++)
    {
        result->uniqueTime[i] = Gt->uniqueTime[nGeTime + i - 1];
        result->survivalProb[i] = Gt->survivalProb[nGeTime + i - 1] / Gt->survivalProb[nGeTime - 1];
    }

    return result;
}

double PseudoRiskTime(
    KMResult *wt,
    double low,
    double high)
{
    if (high < low)
    {
        PRINT_LOCATION();
        printf("Error: high < low when calculating pseudo risk time!\n");
        return (EXIT_FAILURE);
    }

    if (low < 0)
    {
        PRINT_LOCATION();
        printf("Error: lower bound is less than 0 when calculating pseudo risk time!\n");
        return (EXIT_FAILURE);
    }

    double pseudoRiskTime = 0.0;

    if (low > (wt->uniqueTime[wt->nTime - 1] + 1e-8))
    {
        pseudoRiskTime = (high - low) * wt->survivalProb[wt->nTime - 1];
        return pseudoRiskTime;
    }

    for (int i = 1; i < wt->nTime; i++)
    {
        if (wt->uniqueTime[i] < (low + 1e-8))
        {
            continue;
        }
        if ((wt->uniqueTime[i - 1] + 1e-8) > high)
        {
            break;
        }

        double lowTime = fmax(wt->uniqueTime[i - 1], low);
        double highTime = fmin(wt->uniqueTime[i], high);

        pseudoRiskTime += (highTime - lowTime) * wt->survivalProb[i - 1];
    }

    if (wt->uniqueTime[wt->nTime - 1] < high)
    {
        pseudoRiskTime += (high - wt->uniqueTime[wt->nTime - 1]) * wt->survivalProb[wt->nTime - 1];
    }

    return pseudoRiskTime;
}

double PseudoRiskTimeNew(
    KMResult *wt,
    double low,
    double high)
{
    if (high < low)
    {
        PRINT_LOCATION();
        printf("Error: high < low when calculating pseudo risk time!\n");
        return (EXIT_FAILURE);
    }

    if (low < 0)
    {
        PRINT_LOCATION();
        printf("Error: lower bound is less than 0 when calculating pseudo risk time!\n");
        return (EXIT_FAILURE);
    }

    double pseudoRiskTime = 0.0;

    int nsteps = 0;
    double *uniqueTimeNew = (double *)malloc((wt->nTime + 2) * sizeof(double));
    double *survivalProbNew = (double *)malloc((wt->nTime + 2) * sizeof(double));
    uniqueTimeNew[0] = low;
    survivalProbNew[0] = StepFunction(wt->uniqueTime, wt->survivalProb, wt->nTime, low);
    for (int i = 0; i < wt->nTime; i++)
    {
        if (wt->uniqueTime[i] > low && wt->uniqueTime[i] < high)
        {
            nsteps++;
            uniqueTimeNew[nsteps] = wt->uniqueTime[i];
            survivalProbNew[nsteps] = wt->survivalProb[i];
        }
    }
    nsteps++;
    uniqueTimeNew[nsteps] = high;
    survivalProbNew[nsteps] = StepFunction(wt->uniqueTime, wt->survivalProb, wt->nTime, high);

    for (int i = 0; i < nsteps; i++)
    {
        pseudoRiskTime += (uniqueTimeNew[i + 1] - uniqueTimeNew[i]) * survivalProbNew[i];
    }
    free(uniqueTimeNew);
    free(survivalProbNew);
    return pseudoRiskTime;
}
