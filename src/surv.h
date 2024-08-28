
#ifndef SURV_H
#define SURV_H

#include "utils.h"

// Define a struct to hold the result
// *uniqueTime: unique time points where the event happened
// *survivalProb: survival probability at each unique time point
// nTime: number of unique time points
typedef struct _KMResult {
    double *uniqueTime;
    double *survivalProb;
    size_t nTime;
} KMResult;

double StepFunction(double *x, double *y, size_t len, double x0);

// free the memory allocated for KMResult
void FreeKMResult(KMResult *result);


// Kaplan-Meier estimation
//' @param time time of interest
//' @param status status of interest, either 0 or 1, 0 for censored, 1 for terminal event
//' @param nObs number of observations
KMResult *KaplanMeier(double *time, double *status, int nObs);

// wt calculation given
//' @param Gt Kaplan-Meier result w.r.t censoring time
//' @param time time of interest
//' @param status status of interest, either 0 or 1, 0 for censored, 1 for terminal event
KMResult *Wt(KMResult *Gt, double time, double status);

// pseudo risk time calculation
//' @param wt inverse probability of censoring weight w.r.t time of interest
//' @param low lower bound of the pseudo risk time
//' @param high upper bound of the pseudo risk time
double PseudoRiskTime(KMResult *wt, double low, double high);

double PseudoRiskTimeNew(
    KMResult *wt,
    double low,
    double high);

#endif
