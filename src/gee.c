#include "gee.h"

int _GSL_ERROR_FLAG = 0;

int Cgee(
    MATRIX **xin,
    MATRIX **yin,
    MATRIX **idin,
    MATRIX **nin,
    MATRIX **offsetin,
    int *nobs,
    int *p,
    int *parmvec,
    int *M_parm,
    MATRIX **betain,  // keep
    MATRIX **naivvar, // keep
    MATRIX **robvar,  // keep
    MATRIX **Rin,     // keep
    double *S_phi,    // keep
    double *tol,
    int *maxsz,  // keep
    int *S_iter, // keep
    int *silent,
    int *scale_fix,
    int *compatflag)
{
    /* MAIN DECLS */
    if (!(*silent))
    {
        printf("\n nobs = %d", *nobs);
        printf("\n p = %d", *p);
        printf("\n parmvec = %d %d %d", *parmvec, *(parmvec + 1), *(parmvec + 2));
        printf("\nM_parm = %d", *M_parm);
        printf("\nS_phi = %.3f", *S_phi);
        printf("\ntol = %.3f", *tol);
        printf("\nmaxsz = %d", *maxsz);
        printf("\nS_iter = %d", *S_iter);
        printf("\nsilent = %d", *silent);
        printf("\nscale_fix = %d", *scale_fix);
        printf("\ncompatflag = %d", *compatflag);
    }
    MATRIX **X, **Y, **OFFSET;
    MATRIX **N;
    MATRIX *betasave = NULL;
    MATRIX *alpha = NULL;
    MATRIX *R = NULL, *tmpR = NULL, *Ri = NULL, *updateR = NULL;
    MATRIX *One = NULL, *mui = NULL; /*, *Opmui; */
    MATRIX *Ai = NULL, *ei = NULL, *ete = NULL;
    MATRIX *S1 = NULL, *S2 = NULL, *Di = NULL, *this_R = NULL, *S5 = NULL, *S2i = NULL;
    MATRIX *tempmat1 = NULL, *tempmat2 = NULL, *tempmat3 = NULL, *tempmat4 = NULL, *tempmat5 = NULL, *tempmat6 = NULL, *tempmat7 = NULL, *tempmat8 = NULL;
    MATRIX *Aop = NULL, *Dop = NULL, *zi = NULL, *DRop = NULL;
    MATRIX *lag_wts = NULL, *tmpeep = NULL, *scratch = NULL, *wt = NULL;
    double phi, dni, phiLZ;
    int iter, ini, i2, j2, k;
    int alpha_VC_GEE_bandwidth;
    int *onep, one, nclust, i; /*, j; */
    double alpha_scalar;
    double alpha_scalar_LZ, exdiv_LZ;
    int maxni, ni, *maxnip;
    int link, var_mean_rel, corstruct;
    int *maxiter;
    double nnsclust = 0.; /* for counting non-singletons */

    maxiter = (int *)malloc(sizeof(int));

    if (!(*silent))
        printf("@(#) ugee.c 98/01/26 Cgee: GEE C source version chanlib 4.12 \n");

    /* Initialize data */

    alpha = NULL;
    alpha_scalar = 0.;
    alpha_scalar_LZ = 0.;
    exdiv_LZ = 0;
    iter = 0;
    one = 1.;

    *maxiter = *S_iter;
    link = (*parmvec) - 1;
    var_mean_rel = (*(parmvec + 1)) - 1;
    corstruct = (*(parmvec + 2)) - 1;

    // printf("link = %d, var_mean_rel = %d, corstruct = %d\n", link, var_mean_rel, corstruct);

    alpha_VC_GEE_bandwidth = *M_parm + 1;

    onep = &one;

    nclust = VC_GEE_nchanges(idin);
    if (!(*silent))
    {
        printf("xin matrix:\n");
        PrintMatrix(xin);
        printf("yin matrix:\n");
        PrintMatrix(yin);
        printf("idin matrix:\n");
        PrintMatrix(idin);
        printf("offsetin matrix:\n");
        PrintMatrix(offsetin);
        printf("nin matrix:\n");
        PrintMatrix(nin);
        printf("beta in Matrix is\n");
        PrintMatrix(betain);
        printf("R in Matrix is\n");
        PrintMatrix(Rin);
    }
#define set_matrix_array(arrname, nel)                                                           \
    if (!(arrname = (MATRIX **)malloc((unsigned)(nel * sizeof(MATRIX *)))))                      \
    {                                                                                            \
        printf("GEE Error: set_matrix_array (mac): out of memory, requesting %d elements", nel); \
        Free(maxiter);                                                                           \
        return (EXIT_FAILURE);                                                                   \
    }

    set_matrix_array(X, nclust)
        set_matrix_array(Y, nclust)
            set_matrix_array(OFFSET, nclust)
                set_matrix_array(N, nclust)

                    VC_GEE_split(xin, idin, &X);
    VC_GEE_split(yin, idin, &Y);
    VC_GEE_split(offsetin, idin, &OFFSET);
    VC_GEE_split(nin, idin, &N);

    maxni = Y[0]->nrows;
    for (i = 1; i < nclust; i++)
    {
        ni = Y[i]->nrows;
        if (ni > maxni)
            maxni = ni;
    }
    *maxsz = maxni;
    maxnip = &maxni;

    if (corstruct == (int)fixed)
    {
        make_ephemeral(*Rin);
        R = VC_GEE_matcopy(Rin);
    }

    if (!(*silent))
    {
        printf("Cgee will use: ");
        switch (link)
        {
        case VC_GEE_identity:
            printf("VC_GEE_identity link, ");
            break;
        case logarithm:
            printf("log link, ");
            break;
        case logit:
            printf("logit link, ");
            break;
        case reciprocal:
            printf("recip link, ");
            break;
        case probit:
            printf("probit link, ");
            break;
        case cloglog:
            printf("cloglog link, ");
            break;
        default:
            printf("GEE Error: unknown link");
            VC_GEE_destroy_double_matrix(&X, nclust);
            VC_GEE_destroy_double_matrix(&Y, nclust);
            VC_GEE_destroy_double_matrix(&N, nclust);
            VC_GEE_destroy_double_matrix(&OFFSET, nclust);
            VC_GEE_destroy_matrix(&R);
            Free(maxiter);
            return (EXIT_FAILURE);
            break;
        }
        switch (var_mean_rel)
        {
        case Gaussian:
            printf("Gaussian var, ");
            break;
        case Poisson:
            printf("Poisson var, ");
            break;
        case Binomial:
            printf("Binomial var, ");
            break;
        case Gamma:
            printf("Gamma var, ");
            break;
        default:
            printf("GEE Error: Cgee: unknown var_mean_rel. Dies.");
            VC_GEE_destroy_double_matrix(&X, nclust);
            VC_GEE_destroy_double_matrix(&Y, nclust);
            VC_GEE_destroy_double_matrix(&N, nclust);
            VC_GEE_destroy_double_matrix(&OFFSET, nclust);
            VC_GEE_destroy_matrix(&R);
            Free(maxiter);
            return (EXIT_FAILURE);
            break;
        }
        switch (corstruct)
        {
        case independence:
            printf("indep corstr. ");
            break;
        case exchangeable:
            printf("exch corstr. ");
            break;
        case stat_M_dep:
            printf("stat corstr. ");
            break;
        case AR_M:
            printf("AR-M corstr. ");
            break;
        case non_stat_M_dep:
            printf("non-statM corstr. ");
            break;
        case unstructured:
            printf("unstr corstr. ");
            break;
        case fixed:
            printf("fixed corstr. ");
            break;
        default:
            printf("GEE Error: unknown corstr");
            VC_GEE_destroy_double_matrix(&X, nclust);
            VC_GEE_destroy_double_matrix(&Y, nclust);
            VC_GEE_destroy_double_matrix(&N, nclust);
            VC_GEE_destroy_double_matrix(&OFFSET, nclust);
            VC_GEE_destroy_matrix(&R);
            Free(maxiter);
            return (EXIT_FAILURE);
            break;
        }
        printf("\n");
    }

    scratch = VC_GEE_create_matrix(maxni, maxni, EPHEMERAL);
    if (!(*silent))
    {
        printf("scratch matrix:\n");
        PrintMatrix(&scratch);
    }

    updateR = VC_GEE_create_matrix(maxni, maxni, PERMANENT);
    do
    {
        betasave = VC_GEE_matcopy(betain);
        phi = 0.;
        phiLZ = 0.;
        nnsclust = 0.;

        // tmpR = VC_GEE_create_matrix(maxni, maxni, EPHEMERAL);

        switch (corstruct)
        {
        case independence:
            break;
        case exchangeable:
            alpha_scalar = 0.;
            alpha_scalar_LZ = 0.;
            exdiv_LZ = 0.;
            break;
        case stat_M_dep:
        case AR_M:
            alpha = VC_GEE_create_matrix(1, alpha_VC_GEE_bandwidth, EPHEMERAL);
            break;
        case non_stat_M_dep:
        case unstructured:
            alpha = VC_GEE_create_matrix(maxni, maxni, EPHEMERAL);
            break;
        case fixed:
            break;
        default:
            break;
        }

        make_permanent(*betain);
        for (i = 0; i < nclust; i++)
        {
            ni = Y[i]->nrows;
            dni = (double)ni;
            One = VC_GEE_col_1s(ni);
            make_permanent(One);
            tempmat3 = VC_GEE_matmult(&(X[i]), betain);
            tempmat1 = VC_GEE_matadd(&tempmat3, &(OFFSET[i]));
            if (!(*silent))
            {
                printf("****************%d-th cluster****************\n", i);
                printf("X[%d] matrix:\n", i);
                PrintMatrix(&(X[i]));
                printf("OFFSET[%d] matrix:\n", i);
                PrintMatrix(&(OFFSET[i]));
                printf("tempmat1 matrix:\n");
                PrintMatrix(&tempmat1);
            }
            /******************link*******************/
            switch (link)
            {
                double maxfitted;
            case VC_GEE_identity:
                mui = VC_GEE_matcopy(&tempmat1);
                VC_GEE_destroy_matrix(&tempmat1);
                break;
            case logarithm:
                mui = VC_GEE_matexp(&tempmat1);
                break;
            case logit:
                tempmat1 = VC_GEE_matexp(&tempmat1);
                make_permanent(tempmat1);
                tempmat2 = VC_GEE_matadd(&One, &tempmat1);
                tempmat3 = VC_GEE_px1_times_pxq(&(N[i]), &tempmat1);
                mui = VC_GEE_pxq_divby_px1(&tempmat3, &tempmat2);
                make_permanent(mui);
                tempmat3 = VC_GEE_pxq_divby_px1(&mui, &(N[i]));
                maxfitted = VC_GEE_matmax(&tempmat3);
                if ((maxfitted >= .9999 || maxfitted <= .0001) && var_mean_rel == (int)Binomial)
                {
                    printf("GEE Error: Cgee: error: logistic model for probability has fitted value very close to 1.\nestimates diverging; iteration terminated.");
                    Free(maxiter);
                    VC_GEE_destroy_double_matrix(&X, nclust);
                    VC_GEE_destroy_double_matrix(&Y, nclust);
                    VC_GEE_destroy_double_matrix(&N, nclust);
                    VC_GEE_destroy_double_matrix(&OFFSET, nclust);
                    VC_GEE_destroy_matrices(8, &R, &scratch, &updateR, &betasave, &One, &tempmat1, &tempmat2, &mui);
                    return (EXIT_FAILURE);
                }
                break;
            case reciprocal:
                mui = VC_GEE_pxq_divby_px1(&One, &tempmat1);
                break;
                /* probit case added by pj catalano*/ /* OK */
            case probit:
                tempmat1 = VC_GEE_matncdf(&tempmat1);
                mui = VC_GEE_px1_times_pxq(&(N[i]), &tempmat1);
                make_permanent(mui);
                tempmat2 = VC_GEE_pxq_divby_px1(&mui, &(N[i]));
                maxfitted = VC_GEE_matmax(&tempmat2);
                //		if ((maxfitted >= .9999) && var_mean_rel == Binomial)
                if ((maxfitted >= .9999 || maxfitted <= .0001) && var_mean_rel == (int)Binomial)
                {
                    printf("GEE Error: Cgee: estimates diverging; iteration terminated");
                    Free(maxiter);
                    VC_GEE_destroy_double_matrix(&X, nclust);
                    VC_GEE_destroy_double_matrix(&Y, nclust);
                    VC_GEE_destroy_double_matrix(&N, nclust);
                    VC_GEE_destroy_double_matrix(&OFFSET, nclust);
                    VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
                    return (EXIT_FAILURE);
                }
                break;

            /* cloglog case added by jh maindonald */
            case cloglog:
                tempmat1 = VC_GEE_matanticlog(&tempmat1);
                mui = VC_GEE_px1_times_pxq(&(N[i]), &tempmat1);
                make_permanent(mui);
                tempmat3 = VC_GEE_pxq_divby_px1(&mui, &(N[i]));
                maxfitted = VC_GEE_matmax(&tempmat3);
                //		if ((maxfitted >= .999999) && var_mean_rel == (int) Binomial)
                if ((maxfitted >= .9999 || maxfitted <= .0001) && var_mean_rel == (int)Binomial)
                {
                    printf("GEE Error: Cgee: error: cloglog model for probability has fit ted value very close to 1.\nestimates diverging; iteration terminated.");
                    Free(maxiter);
                    VC_GEE_destroy_double_matrix(&X, nclust);
                    VC_GEE_destroy_double_matrix(&Y, nclust);
                    VC_GEE_destroy_double_matrix(&N, nclust);
                    VC_GEE_destroy_double_matrix(&OFFSET, nclust);
                    VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
                    return (EXIT_FAILURE);
                }
                break;

            default:
                printf("GEE Error: Cgee: unknown link. Dies.");
                Free(maxiter);
                VC_GEE_destroy_double_matrix(&X, nclust);
                VC_GEE_destroy_double_matrix(&Y, nclust);
                VC_GEE_destroy_double_matrix(&N, nclust);
                VC_GEE_destroy_double_matrix(&OFFSET, nclust);
                VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
                return (EXIT_FAILURE);
                break;
            }
            make_permanent(mui);
            ei = VC_GEE_matsub(&(Y[i]), &mui);
            if (!(*silent))
            {
                printf("mui matrix:\n");
                PrintMatrix(&mui);
                printf("ei matrix:\n");
                PrintMatrix(&ei);
            }
            switch (var_mean_rel)
            {
            case Gaussian:
                /*	Ai = VC_GEE_ident(ni); */
                break;
            case Poisson:
                Ai = VC_GEE_form_diag(&mui);
                break;
            case Binomial:
                tempmat3 = VC_GEE_pxq_divby_px1(&mui, &(N[i]));
                tempmat4 = VC_GEE_matsub(&One, &tempmat3);
                tempmat5 = VC_GEE_px1_times_pxq(&mui, &tempmat4);
                Ai = VC_GEE_form_diag(&tempmat5);
                break;
            case Gamma:
                tempmat3 = VC_GEE_px1_times_pxq(&mui, &mui);
                Ai = VC_GEE_form_diag(&tempmat3);
                break;
            default:
                printf("GEE Error: Cgee: unknown var_mean_rel. Dies.");
                Free(maxiter);
                VC_GEE_destroy_double_matrix(&X, nclust);
                VC_GEE_destroy_double_matrix(&Y, nclust);
                VC_GEE_destroy_double_matrix(&N, nclust);
                VC_GEE_destroy_double_matrix(&OFFSET, nclust);
                VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
                return (EXIT_FAILURE);
                break;
            }

            VC_GEE_destroy_matrix(&mui);

            if (var_mean_rel != Gaussian)
            {
                tempmat3 = VC_GEE_diag_as_vec(&Ai);
                tempmat4 = VC_GEE_matsqrt(&tempmat3);
                tempmat5 = VC_GEE_mat1over(&tempmat4);
                ei = VC_GEE_px1_times_pxq(&tempmat5, &ei);
            }
            make_permanent(ei);
            tempmat3 = VC_GEE_transp(&ei);
            ete = VC_GEE_matmult(&tempmat3, &ei);

            phi += ete->data[0][0] / dni;
            phiLZ += ete->data[0][0];
            if (dni > 1)
                nnsclust += 1.0;

            if (corstruct == (int)stat_M_dep || corstruct == (int)AR_M)
            {
                if (dni < (double)alpha_VC_GEE_bandwidth)
                {
                    printf("GEE Error: cgee: M-dependence, M=%d, but clustsize=%d\nfatal error for this model", (int)*M_parm, (int)dni);
                    Free(maxiter);
                    VC_GEE_destroy_double_matrix(&X, nclust);
                    VC_GEE_destroy_double_matrix(&Y, nclust);
                    VC_GEE_destroy_double_matrix(&N, nclust);
                    VC_GEE_destroy_double_matrix(&OFFSET, nclust);
                    VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
                    return (EXIT_FAILURE);
                }
                lag_wts = VC_GEE_create_matrix(1, alpha_VC_GEE_bandwidth, EPHEMERAL);
                lag_wts->data[0][0] = (double)1.;

                for (ini = 1; ini < alpha_VC_GEE_bandwidth; ini++)
                {
                    lag_wts->data[0][ini] = dni / (dni - (double)ini);
                }
            }

            switch (corstruct)
            {
            case independence:
                break;
            case exchangeable:
                if (ni > 1)
                    tempmat3 = VC_GEE_transp(&ei);
                tempmat4 = VC_GEE_matmult(&ei, &tempmat3);
                alpha_scalar +=
                    1. / (dni * (dni - 1.)) * (VC_GEE_elsum(&tempmat4) - ete->data[0][0]);
                if (ni > 1)
                    tempmat3 = VC_GEE_transp(&ei);
                tempmat4 = VC_GEE_matmult(&ei, &tempmat3);
                alpha_scalar_LZ +=
                    (VC_GEE_elsum(&tempmat4) - ete->data[0][0]);
                if (ni > 1)
                    exdiv_LZ += dni * (dni - 1); /* suppress .5 because num
                                  is redundant */
                break;
            case stat_M_dep:
            case AR_M:
                tempmat3 = VC_GEE_covlag(&ei, alpha_VC_GEE_bandwidth, 0);
                tempmat4 = VC_GEE_transp(&tempmat3);
                tempmat5 = VC_GEE_transp(&lag_wts);
                tempmat6 = VC_GEE_px1_times_pxq(&tempmat5, &tempmat4);
                tempmat7 = VC_GEE_transp(&tempmat6);
                alpha = VC_GEE_matadd(&alpha, &tempmat7);
                // alpha = VC_GEE_matadd(alpha, VC_GEE_transp(VC_GEE_px1_times_pxq(VC_GEE_transp(lag_wts), VC_GEE_transp(VC_GEE_covlag(ei, alpha_VC_GEE_bandwidth, 0)))));
                break;
            case non_stat_M_dep:
            case unstructured:
                tempmat3 = VC_GEE_transp(&ei);
                tmpeep = VC_GEE_matmult(&ei, &tempmat3);
                VC_GEE_plug(&tmpeep, &scratch, 0, 0);
                make_permanent(scratch);
                alpha = VC_GEE_matadd(&alpha, &scratch);
                make_ephemeral(scratch);
                scratch = VC_GEE_scalar_times_matrix(0, &scratch);
                break;
            case fixed:
                break;
            default:
                printf("GEE Error: corstruct not implemented.");
                Free(maxiter);
                VC_GEE_destroy_double_matrix(&X, nclust);
                VC_GEE_destroy_double_matrix(&Y, nclust);
                VC_GEE_destroy_double_matrix(&N, nclust);
                VC_GEE_destroy_double_matrix(&OFFSET, nclust);
                VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
                return (EXIT_FAILURE);
                break;
            }
            VC_GEE_destroy_matrices(3, &ei, &One, &ete);
            // VC_GEE_destroy_matrix(&ei);
            // VC_GEE_destroy_matrix(&One);
            // VC_GEE_destroy_matrix(&ete);
        }

        if (alpha != NULL)
            alpha = VC_GEE_scalar_times_matrix(
                (double)nclust / (phi * nnsclust), &alpha);

        alpha_scalar /= phi;
        alpha_scalar *= (double)nclust / nnsclust;

        phi /= (double)nclust;
        phiLZ /= ((double)(*nobs - *p));

        alpha_scalar_LZ /= (phiLZ * (exdiv_LZ - 2. * (double)*p));
        if (*compatflag == 0)
        { /* abandon compatibility with early macro */
            alpha_scalar = alpha_scalar_LZ;
            phi = phiLZ;
        }

        switch (corstruct)
        {
        case independence: /* this is ridiculous */
            R = VC_GEE_ident(maxni);
            break;
        case exchangeable:
            tempmat3 = VC_GEE_col_1s(maxni);
            make_permanent(tempmat3);
            tempmat4 = VC_GEE_transp(&tempmat3);
            tempmat5 = VC_GEE_matmult(&tempmat3, &tempmat4);
            tempmat6 = VC_GEE_scalar_times_matrix(alpha_scalar, &tempmat5);
            VC_GEE_destroy_matrix(&tempmat3);

            tempmat3 = VC_GEE_ident(maxni);
            tempmat4 = VC_GEE_scalar_times_matrix((double)1. - alpha_scalar, &tempmat3);
            R = VC_GEE_matadd(&tempmat6, &tempmat4);
            break;
        case stat_M_dep:
            tmpR = VC_GEE_create_matrix(1, maxni, EPHEMERAL);
            VC_GEE_plug(&alpha, &tmpR, 0, 0);
            tmpR->data[0][0] = (double)1.;
            R = VC_GEE_toeplitz(&tmpR);
            break;
        case AR_M:
            tmpR = VC_GEE_create_matrix(1, maxni, EPHEMERAL);
            alpha->data[0][0] = (double)1.;
            make_permanent(alpha);
            VC_GEE_plug(&alpha, &tmpR, 0, 0);
            tempmat1 = VC_GEE_extract_cols(&alpha, 0, alpha_VC_GEE_bandwidth - 2);
            wt = VC_GEE_toeplitz(&tempmat1);
            tempmat1 = VC_GEE_extract_cols(&alpha, 1, alpha_VC_GEE_bandwidth - 1);
            VC_GEE_destroy_matrix(&alpha);
            tempmat3 = matrix_inverse(&wt);
            if (tempmat3 == NULL)
            {
                if (_GSL_ERROR_FLAG)
                    printf("Warning: matrix inverse fail in %s:%d, nobs = %d\n", __FILE__, __LINE__ - 4, *nobs);
                Free(maxiter);
                VC_GEE_destroy_double_matrix(&X, nclust);
                VC_GEE_destroy_double_matrix(&Y, nclust);
                VC_GEE_destroy_double_matrix(&N, nclust);
                VC_GEE_destroy_double_matrix(&OFFSET, nclust);
                VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
                return (INV_FAILURE);
            }
            alpha = VC_GEE_matmult(&tempmat1, &tempmat3);
            for (i2 = alpha_VC_GEE_bandwidth; i2 < maxni; i2++)
            {
                for (j2 = 0; j2 < alpha_VC_GEE_bandwidth - 1; j2++)
                {
                    tmpR->data[0][i2] += alpha->data[0][j2] * tmpR->data[0][(i2 - j2) - 1];
                }
            }
            R = VC_GEE_toeplitz(&tmpR);
            VC_GEE_destroy_matrix(&alpha);
            break;
        case non_stat_M_dep:
            R = VC_GEE_band(&alpha, alpha_VC_GEE_bandwidth);
            for (k = 0; k < R->ncols; k++)
            {
                R->data[k][k] = (double)1.;
            }
            break;
        case unstructured:
            R = VC_GEE_matcopy(&alpha);
            for (k = 0; k < R->ncols; k++)
            {
                R->data[k][k] = (double)1.;
            }
            VC_GEE_destroy_matrix(&alpha);
            break;
        case fixed:
            break;
        default:
            printf("GEE Error: corstruct not implemented.");
            Free(maxiter);
            VC_GEE_destroy_double_matrix(&X, nclust);
            VC_GEE_destroy_double_matrix(&Y, nclust);
            VC_GEE_destroy_double_matrix(&N, nclust);
            VC_GEE_destroy_double_matrix(&OFFSET, nclust);
            VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
            return (EXIT_FAILURE);
            break;
        }
        make_permanent(R);
        if (!(*silent))
        {
            printf("R matrix:\n");
            PrintMatrix(&R);
        }

        S1 = VC_GEE_create_matrix((int)*p, 1, EPHEMERAL);
        S2 = VC_GEE_create_matrix((int)*p, (int)*p, EPHEMERAL);

        for (i = 0; i < nclust; i++)
        {
            ni = Y[i]->nrows;
            dni = (double)ni;
            One = VC_GEE_col_1s(ni);
            make_permanent(One);

            tempmat3 = VC_GEE_matmult(&(X[i]), betain);
            tempmat1 = VC_GEE_matadd(&tempmat3, &(OFFSET[i]));
            if (!(*silent))
            {
                printf("************%d-th cluster second loop*************\n", i);
                printf("tempmat1 matrix:\n");
                PrintMatrix(&tempmat1);
            }
            switch (link)
            {
            case VC_GEE_identity:
                Di = VC_GEE_matcopy(&(X[i]));
                mui = VC_GEE_matcopy(&tempmat1);
                VC_GEE_destroy_matrix(&tempmat1);
                break;
            case logarithm:
                tempmat1 = VC_GEE_matexp(&tempmat1);
                mui = VC_GEE_matcopy(&tempmat1);
                Di = VC_GEE_px1_times_pxq(&tempmat1, &(X[i]));
                break;
            case logit:
                tempmat1 = VC_GEE_matexp(&tempmat1);
                make_permanent(tempmat1);
                tempmat2 = VC_GEE_matadd(&One, &tempmat1);
                tempmat3 = VC_GEE_px1_times_pxq(&(N[i]), &tempmat1);
                tempmat2 = VC_GEE_pxq_divby_px1(&tempmat3, &tempmat2);
                make_permanent(tempmat2);
                VC_GEE_destroy_matrix(&tempmat1);
                tempmat3 = VC_GEE_pxq_divby_px1(&tempmat2, &(N[i]));
                tempmat1 = VC_GEE_matsub(&One, &tempmat3);
                tempmat1 = VC_GEE_px1_times_pxq(&tempmat1, &tempmat2);
                Di = VC_GEE_px1_times_pxq(&tempmat1, &(X[i]));
                mui = VC_GEE_matcopy(&tempmat2);
                VC_GEE_destroy_matrix(&tempmat2);
                break;
            case reciprocal:
                tempmat1 = VC_GEE_pxq_divby_px1(&One, &tempmat1);
                mui = VC_GEE_matcopy(&tempmat1);
                tempmat2 = VC_GEE_matcopy(&tempmat1);
                tempmat1 = VC_GEE_px1_times_pxq(&tempmat1, &tempmat2);
                tempmat1 = VC_GEE_px1_times_pxq(&tempmat1, &(X[i]));
                Di = VC_GEE_scalar_times_matrix(-1., &tempmat1);
                break;
            /* probit case added by pj catalano */
            case probit:
                tempmat2 = VC_GEE_matcopy(&tempmat1);
                tempmat2 = VC_GEE_matncdf(&tempmat2);
                mui = VC_GEE_px1_times_pxq(&(N[i]), &tempmat2);
                tempmat1 = VC_GEE_matnpdf(&tempmat1);
                tempmat2 = VC_GEE_px1_times_pxq(&(N[i]), &tempmat1);
                Di = VC_GEE_px1_times_pxq(&tempmat2, &(X[i]));
                break;
                /* cloglog case added by jh maindonald */
            case cloglog:
                tempmat2 = VC_GEE_matcopy(&tempmat1);
                tempmat1 = VC_GEE_matanticlog(&tempmat1);
                make_permanent(tempmat1);
                mui = VC_GEE_px1_times_pxq(&(N[i]), &(tempmat1));
                tempmat2 = VC_GEE_matexp(&(tempmat2));
                tempmat2 = VC_GEE_px1_times_pxq(&(N[i]), &tempmat2);
                tempmat3 = VC_GEE_matsub(&One, &tempmat1);
                tempmat1 = VC_GEE_px1_times_pxq(&tempmat2, &tempmat3);
                Di = VC_GEE_px1_times_pxq(&tempmat1, &(X[i]));
                VC_GEE_destroy_matrix(&tempmat1);
                break;

            default:
                printf("GEE Error: Cgee: unknown link. Dies.");
                Free(maxiter);
                VC_GEE_destroy_double_matrix(&X, nclust);
                VC_GEE_destroy_double_matrix(&Y, nclust);
                VC_GEE_destroy_double_matrix(&N, nclust);
                VC_GEE_destroy_double_matrix(&OFFSET, nclust);
                VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
                return (EXIT_FAILURE);
                break;
            }
            make_permanent(mui);
            make_permanent(Di);

            tempmat3 = VC_GEE_matsub(&(Y[i]), &mui);
            tempmat4 = VC_GEE_matmult(&Di, betain);
            zi = VC_GEE_matadd(&tempmat4, &tempmat3);

            switch (var_mean_rel)
            {
            case Gaussian:
                /* Ai = VC_GEE_ident(ni); */
                break;
            case Poisson:
                Ai = VC_GEE_form_diag(&mui);
                break;
            case Binomial:
                tempmat3 = VC_GEE_pxq_divby_px1(&mui, &(N[i]));
                tempmat4 = VC_GEE_matsub(&One, &tempmat3);
                tempmat5 = VC_GEE_px1_times_pxq(&mui, &tempmat4);
                Ai = VC_GEE_form_diag(&tempmat5);
                // Ai = VC_GEE_form_diag(VC_GEE_px1_times_pxq(
                //     mui, VC_GEE_matsub(One, VC_GEE_pxq_divby_px1(mui, N[i]))));
                break;
            case Gamma:
                tempmat3 = VC_GEE_px1_times_pxq(&mui, &mui);
                Ai = VC_GEE_form_diag(&tempmat3);
                break;
            default:
                printf("GEE Error: Cgee: unknown var_mean_rel. Dies.");
                Free(maxiter);
                VC_GEE_destroy_double_matrix(&X, nclust);
                VC_GEE_destroy_double_matrix(&Y, nclust);
                VC_GEE_destroy_double_matrix(&N, nclust);
                VC_GEE_destroy_double_matrix(&OFFSET, nclust);
                VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
                return (EXIT_FAILURE);
                break;
            }
            VC_GEE_destroy_matrix(&mui);

            if (var_mean_rel != Gaussian) /* else Ai is VC_GEE_identity */
            {
                tempmat3 = VC_GEE_diag_as_vec(&Ai);
                tempmat4 = VC_GEE_matsqrt(&tempmat3);
                Aop = VC_GEE_mat1over(&tempmat4);
                // Aop = VC_GEE_mat1over(VC_GEE_matsqrt(VC_GEE_diag_as_vec(Ai)));
            }
            if (var_mean_rel != Gaussian) /* else Ai is VC_GEE_identity */
                make_permanent(Aop);

            make_ephemeral(Di);
            if (var_mean_rel != Gaussian)
            {
                // printf("I'm here\n");
                Di = VC_GEE_px1_times_pxq(&Aop, &Di);
            }
            make_permanent(Di);
            if (var_mean_rel != Gaussian)
                zi = VC_GEE_px1_times_pxq(&Aop, &zi);
            if (!(*silent))
            {
                printf("Di matrix:\n");
                PrintMatrix(&Di);
                printf("zi matrix:\n");
                PrintMatrix(&zi);
            }
            if (corstruct != independence)
            {
                this_R = VC_GEE_corner(&R, ni, ni);
                Ri = matrix_inverse(&this_R);
                if (Ri == NULL)
                {
                    if (_GSL_ERROR_FLAG)
                        printf("Warning: matrix inverse fail in %s:%d, nobs = %d\n", __FILE__, __LINE__ - 4, *nobs);
                    Free(maxiter);
                    VC_GEE_destroy_double_matrix(&X, nclust);
                    VC_GEE_destroy_double_matrix(&Y, nclust);
                    VC_GEE_destroy_double_matrix(&N, nclust);
                    VC_GEE_destroy_double_matrix(&OFFSET, nclust);
                    VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
                    return (INV_FAILURE);
                }
                tempmat3 = VC_GEE_transp(&Di);
                Dop = VC_GEE_matmult(&tempmat3, &Ri);
            }
            else
            {
                Dop = VC_GEE_transp(&Di);
            }

            make_permanent(Dop);
            if (!(*silent))
            {
                printf("Dop matrix:\n");
                PrintMatrix(&Dop);
            }
            tempmat3 = VC_GEE_matmult(&Dop, &zi);
            tempmat4 = VC_GEE_matmult(&Dop, &Di);
            S1 = VC_GEE_matadd(&S1, &tempmat3);
            S2 = VC_GEE_matadd(&S2, &tempmat4);

            VC_GEE_destroy_matrix(&Dop);
            if (var_mean_rel != Gaussian)
                VC_GEE_destroy_matrix(&Aop);
            VC_GEE_destroy_matrix(&One);
            VC_GEE_destroy_matrix(&Di); // permanent
        }

        VC_GEE_destroy_matrix(betain);
        tempmat3 = matrix_inverse(&S2);
        if (tempmat3 == NULL)
        {
            if (_GSL_ERROR_FLAG)
                printf("Warning: matrix inverse fail in %s:%d, nobs = %d\n", __FILE__, __LINE__ - 4, *nobs);
            Free(maxiter);
            VC_GEE_destroy_double_matrix(&X, nclust);
            VC_GEE_destroy_double_matrix(&Y, nclust);
            VC_GEE_destroy_double_matrix(&N, nclust);
            VC_GEE_destroy_double_matrix(&OFFSET, nclust);
            VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
            return (INV_FAILURE);
        }
        *betain = VC_GEE_matmult(&tempmat3, &S1);
        make_permanent(*betain);

        if (!(*silent))
        {
            printf("current parameter estimates:\n");
            VC_GEE_matdump(betain);
            printf("********************completed iteration %d****************\n", iter);
        }
        iter++;

        for (int i = 0; i < maxni; i++)
        {
            for (int j = 0; j < maxni; j++)
            {
                updateR->data[i][j] = R->data[i][j];
            }
        }
        if (corstruct != fixed)
        {
            VC_GEE_destroy_matrix(&R);
        }
        tempmat3 = VC_GEE_col_1s((int)*p);
        tempmat4 = VC_GEE_pxq_divby_px1(&betasave, betain);
        tempmat5 = VC_GEE_matsub(&tempmat4, &tempmat3);
        tempmat6 = VC_GEE_matabs(&tempmat5);

    } while ((VC_GEE_matmax(&tempmat6) > *tol) && (iter < *maxiter));

    // if (corstruct == fixed)
    // {
    //     VC_GEE_destroy_matrix(&R);
    // }

    if (iter >= *maxiter)
    {
        if (!(*silent))
        {
            printf("Warning: Maximum number of iterations consumed\n");
            printf("Warning: Convergence not achieved; results suspect\n");
        }
        Free(maxiter);
        VC_GEE_destroy_double_matrix(&X, nclust);
        VC_GEE_destroy_double_matrix(&Y, nclust);
        VC_GEE_destroy_double_matrix(&N, nclust);
        VC_GEE_destroy_double_matrix(&OFFSET, nclust);
        VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
        return ((int)MAXITER_EXCEEDED);
    }

    S2 = VC_GEE_create_matrix((int)*p, (int)*p, EPHEMERAL);
    S5 = VC_GEE_create_matrix((int)*p, (int)*p, EPHEMERAL);
    phi = 0.;
    phiLZ = 0.;

    for (i = 0; i < nclust; i++)
    {
        ni = Y[i]->nrows;
        dni = (double)ni;
        One = VC_GEE_col_1s(ni);
        make_permanent(One);
        tempmat3 = VC_GEE_matmult(&X[i], betain);
        tempmat1 = VC_GEE_matadd(&tempmat3, &OFFSET[i]);
        if (!(*silent))
        {
            printf("************%d-th cluster third loop*************\n", i);
            printf("tempmat1 matrix:\n");
            PrintMatrix(&tempmat1);
        }
        switch (link)
        {
        case VC_GEE_identity:
            Di = VC_GEE_matcopy(&X[i]);
            mui = VC_GEE_matcopy(&tempmat1);
            VC_GEE_destroy_matrix(&tempmat1);
            break;
        case logarithm:
            tempmat1 = VC_GEE_matexp(&tempmat1);
            mui = VC_GEE_matcopy(&tempmat1);
            Di = VC_GEE_px1_times_pxq(&tempmat1, &X[i]);
            break;
        case logit:
            tempmat1 = VC_GEE_matexp(&tempmat1);
            make_permanent(tempmat1);
            tempmat2 = VC_GEE_matadd(&One, &tempmat1);
            tempmat3 = VC_GEE_px1_times_pxq(&N[i], &tempmat1);
            tempmat2 = VC_GEE_pxq_divby_px1(&tempmat3, &tempmat2);
            make_permanent(tempmat2);
            VC_GEE_destroy_matrix(&tempmat1);
            tempmat3 = VC_GEE_pxq_divby_px1(&tempmat2, &N[i]);
            tempmat1 = VC_GEE_matsub(&One, &tempmat3);
            tempmat1 = VC_GEE_px1_times_pxq(&tempmat1, &tempmat2);
            Di = VC_GEE_px1_times_pxq(&tempmat1, &X[i]);
            mui = VC_GEE_matcopy(&tempmat2);
            VC_GEE_destroy_matrix(&tempmat2);
            break;
        case reciprocal:
            tempmat1 = VC_GEE_pxq_divby_px1(&One, &tempmat1);
            mui = VC_GEE_matcopy(&tempmat1);
            tempmat2 = VC_GEE_matcopy(&tempmat1);
            tempmat1 = VC_GEE_px1_times_pxq(&tempmat1, &tempmat2);
            tempmat1 = VC_GEE_px1_times_pxq(&tempmat1, &X[i]);
            Di = VC_GEE_scalar_times_matrix(-1., &tempmat1);
            break;
            /* probit case added by pj catalano*/
        case probit:
            tempmat2 = VC_GEE_matcopy(&tempmat1);
            tempmat2 = VC_GEE_matncdf(&tempmat2);
            mui = VC_GEE_px1_times_pxq(&N[i], &tempmat2);
            tempmat1 = VC_GEE_matnpdf(&tempmat1);
            tempmat2 = VC_GEE_px1_times_pxq(&N[i], &tempmat1);
            Di = VC_GEE_px1_times_pxq(&tempmat2, &X[i]);
            break;
            /* maind */
        case cloglog:
            tempmat2 = VC_GEE_matcopy(&tempmat1);
            tempmat1 = VC_GEE_matanticlog(&tempmat1);
            make_permanent(tempmat1);
            mui = VC_GEE_px1_times_pxq(&N[i], &tempmat1);
            tempmat2 = VC_GEE_matexp(&tempmat2);
            tempmat2 = VC_GEE_px1_times_pxq(&N[i], &tempmat2);
            tempmat3 = VC_GEE_matsub(&One, &tempmat1);
            tempmat1 = VC_GEE_px1_times_pxq(&tempmat2, &tempmat3);
            Di = VC_GEE_px1_times_pxq(&tempmat1, &X[i]);
            VC_GEE_destroy_matrix(&tempmat1);
            break;

        default:
            printf("GEE Error: Cgee: unknown link. Dies.");
            Free(maxiter);
            VC_GEE_destroy_double_matrix(&X, nclust);
            VC_GEE_destroy_double_matrix(&Y, nclust);
            VC_GEE_destroy_double_matrix(&N, nclust);
            VC_GEE_destroy_double_matrix(&OFFSET, nclust);
            VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
            return (EXIT_FAILURE);
            break;
        }
        make_permanent(mui);

        ei = VC_GEE_matsub(&Y[i], &mui);
        if (!(*silent))
        {
            printf("ei matrix:\n");
            PrintMatrix(&ei);
        }
        switch (var_mean_rel)
        {
        case Gaussian:
            /* Ai = VC_GEE_ident(ni); */
            break;
        case Poisson:
            Ai = VC_GEE_form_diag(&mui);
            break;
        case Binomial:
            tempmat3 = VC_GEE_pxq_divby_px1(&mui, &N[i]);
            tempmat4 = VC_GEE_matsub(&One, &tempmat3);
            tempmat5 = VC_GEE_px1_times_pxq(&mui, &tempmat4);
            Ai = VC_GEE_form_diag(&tempmat5);
            break;
        case Gamma:
            tempmat3 = VC_GEE_px1_times_pxq(&mui, &mui);
            Ai = VC_GEE_form_diag(&tempmat3);
            break;
        default:
            printf("GEE Error: Cgee: unknown var_mean_rel. Dies.\n");
            Free(maxiter);
            VC_GEE_destroy_double_matrix(&X, nclust);
            VC_GEE_destroy_double_matrix(&Y, nclust);
            VC_GEE_destroy_double_matrix(&N, nclust);
            VC_GEE_destroy_double_matrix(&OFFSET, nclust);
            VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
            return (EXIT_FAILURE);
            break;
        }
        VC_GEE_destroy_matrices(2, &mui, &One);

        if (var_mean_rel != Gaussian)
        {
            tempmat3 = VC_GEE_diag_as_vec(&Ai);
            tempmat4 = VC_GEE_matsqrt(&tempmat3);
            Aop = VC_GEE_mat1over(&tempmat4);
            make_permanent(Aop);
            Di = VC_GEE_px1_times_pxq(&Aop, &Di);
        }
        make_permanent(Di);

        if (var_mean_rel != Gaussian)
            ei = VC_GEE_px1_times_pxq(&Aop, &ei);
        make_permanent(ei);

        tempmat3 = VC_GEE_transp(&ei);
        ete = VC_GEE_matmult(&tempmat3, &ei);
        phi += ete->data[0][0] / dni;
        phiLZ += ete->data[0][0];

        if (corstruct != independence)
        {
            tempmat3 = VC_GEE_corner(&updateR, ni, ni);
            tempmat4 = matrix_inverse(&tempmat3);
            if (tempmat4 == NULL)
            {
                if (_GSL_ERROR_FLAG)
                    printf("Warning: matrix inverse fail in %s:%d, nobs = %d\n", __FILE__, __LINE__ - 4, *nobs);
                Free(maxiter);
                VC_GEE_destroy_double_matrix(&X, nclust);
                VC_GEE_destroy_double_matrix(&Y, nclust);
                VC_GEE_destroy_double_matrix(&N, nclust);
                VC_GEE_destroy_double_matrix(&OFFSET, nclust);
                VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
                return (INV_FAILURE);
            }
            tempmat5 = VC_GEE_transp(&Di);
            DRop = VC_GEE_matmult(&tempmat5, &tempmat4);
        }
        else
            DRop = VC_GEE_transp(&Di);

        make_permanent(DRop);

        tempmat3 = VC_GEE_matmult(&DRop, &Di);
        S2 = VC_GEE_matadd(&S2, &tempmat3);
        tempmat3 = VC_GEE_transp(&DRop);
        tempmat4 = VC_GEE_transp(&ei);
        tempmat5 = VC_GEE_matmult(&tempmat4, &tempmat3);
        tempmat6 = VC_GEE_matmult(&DRop, &ei);
        tempmat7 = VC_GEE_matmult(&tempmat6, &tempmat5);
        S5 = VC_GEE_matadd(&S5, &tempmat7);

        if (var_mean_rel != Gaussian)
            VC_GEE_destroy_matrix(&Aop);

        VC_GEE_destroy_matrices(4, &ete, &DRop, &Di, &ei);
    }

    phi /= (double)nclust;
    phiLZ /= (double)(*nobs - *p);
    if (*compatflag == 0)
    { /* abandon compatibility with early macro */
        phi = phiLZ;
    }

    S2i = matrix_inverse(&S2);
    if (S2i == NULL)
    {
        if (_GSL_ERROR_FLAG)
            printf("Warning: matrix inverse fail in %s:%d, nobs = %d\n", __FILE__, __LINE__ - 4, *nobs);
        Free(maxiter);
        VC_GEE_destroy_double_matrix(&X, nclust);
        VC_GEE_destroy_double_matrix(&Y, nclust);
        VC_GEE_destroy_double_matrix(&N, nclust);
        VC_GEE_destroy_double_matrix(&OFFSET, nclust);
        VC_GEE_destroy_matrices(25, &betasave, &R, &tmpR, &Ri, &updateR, &One, &mui, &Ai, &ei, &ete, &S1, &S2, &Di, &this_R, &S5, &S2i, &tempmat1, &tempmat2, &Aop, &Dop, &zi, &DRop, &lag_wts, &tmpeep, &scratch, &wt);
        return (INV_FAILURE);
    }
    make_permanent(S2i);
    if (*scale_fix)
        phi = *S_phi;
    if (*scale_fix && var_mean_rel == Gaussian)
        printf("Warning: Scale parameter fixed at %.3f with Gaussian variance function",
               *S_phi);
    *naivvar = VC_GEE_scalar_times_matrix(phi, &S2i);
    tempmat3 = VC_GEE_matmult(&S2i, &S5);
    *robvar = VC_GEE_matmult(&tempmat3, &S2i);

    *S_phi = phi;
    *S_iter = iter;

    make_ephemeral(updateR);
    VC_GEE_destroy_matrix(Rin);
    *Rin = VC_GEE_matcopy(&updateR);
    VC_GEE_destroy_matrix(&updateR);

    Free(maxiter);
    VC_GEE_destroy_double_matrix(&X, nclust);
    VC_GEE_destroy_double_matrix(&Y, nclust);
    VC_GEE_destroy_double_matrix(&N, nclust);
    VC_GEE_destroy_double_matrix(&OFFSET, nclust);
    VC_GEE_destroy_matrices(3, &R, &S2i, &scratch);

    return (EXIT_SUCCESS);
}

/* gee support @(#) chanmat.c 4.12 98/01/26 */
/* @(#) chanmat.nw 1.3 94/03/09 */

MATRIX *VC_GEE_create_matrix(int nrows, int ncols, int permanence)
{
    MATRIX *tmp = (MATRIX *)malloc(sizeof(struct matrix));

    if (tmp == NULL)
    {
        printf("GEE Error: VC_GEE_create_matrix: malloc failed %lu",
               sizeof(struct matrix));
        return (NULL);
    }

    tmp->nrows = nrows;
    tmp->ncols = ncols;
    tmp->permanence = permanence;

    tmp->data = (double **)malloc(nrows * sizeof(double *));
    for (int i = 0; i < nrows; i++)
    {
        tmp->data[i] = (double *)calloc(ncols, sizeof(double));
        if (tmp->data[i] == NULL)
        {
            printf("GEE Error: VC_GEE_create_matrix: malloc failed, nrows=%d ncols=%d",
                   nrows, ncols);
            return (NULL);
        }
    }

    return tmp;
}

MATRIX *VC_GEE_transp(MATRIX **mat)
{
    int i, j;

    MATRIX *tmp = VC_GEE_create_matrix((*mat)->ncols, (*mat)->nrows, EPHEMERAL);
    for (i = 0; i < (*mat)->nrows; i++)
    {
        for (j = 0; j < (*mat)->ncols; j++)
        {
            tmp->data[j][i] = (*mat)->data[i][j];
        }
    }
    free_if_ephemeral(mat);
    return tmp;
}

MATRIX *VC_GEE_corner(MATRIX **mat, int nr, int nc)
{
    if ((nr > (*mat)->nrows) || (nc > (*mat)->ncols))
    {
        printf("GEE Error: VC_GEE_corner: request not a submatrix.\n");
        return (NULL);
    }
    int i, j;
    MATRIX *tmp = VC_GEE_create_matrix(nr, nc, EPHEMERAL);
    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            tmp->data[i][j] = (*mat)->data[i][j];
        }
    }
    free_if_ephemeral(mat);
    return tmp;
}

MATRIX *VC_GEE_extract_rows(MATRIX **Source, int VC_GEE_start, int end)
/* purely zero-based */
{
    if (VC_GEE_start < 0 || end < 0 || VC_GEE_start > (*Source)->nrows || end > (*Source)->nrows || VC_GEE_start > end)
    {
        printf("GEE Error: VC_GEE_extract_rows: invalid start or end");
        return (NULL);
    }
    MATRIX *temp;
    int rows_to_get, i, j;

    rows_to_get = end - VC_GEE_start + 1;

    temp = VC_GEE_create_matrix(rows_to_get, (*Source)->ncols, EPHEMERAL);

    for (i = 0; i < rows_to_get; i++)
    {
        for (j = 0; j < (*Source)->ncols; j++)
        {
            temp->data[i][j] = (*Source)->data[VC_GEE_start + i][j];
        }
    }
    /* DOES NOT CLEAN */
    return temp;
}

MATRIX *VC_GEE_extract_cols(MATRIX **x, int VC_GEE_start, int end)
{
    if (VC_GEE_start < 0 || end < 0 || VC_GEE_start > (*x)->ncols || end > (*x)->ncols || VC_GEE_start > end)
    {
        printf("GEE Error: VC_GEE_extract_cols: invalid start or end");
        return (NULL);
    }
    int cols_to_get = end - VC_GEE_start + 1;
    MATRIX *tmp = VC_GEE_create_matrix((*x)->nrows, cols_to_get, EPHEMERAL);
    int i, j;
    for (i = 0; i < (*x)->nrows; i++)
    {
        for (j = 0; j < cols_to_get; j++)
        {
            tmp->data[i][j] = (*x)->data[i][j + VC_GEE_start];
        }
    }
    free_if_ephemeral(x);
    return tmp;
}

MATRIX *VC_GEE_matcopy(MATRIX **inmat)
{
    int i, j;
    MATRIX *outmat;

    outmat = VC_GEE_create_matrix((*inmat)->nrows, (*inmat)->ncols, EPHEMERAL);
    for (i = 0; i < (*inmat)->nrows; i++)
    {
        for (j = 0; j < (*inmat)->ncols; j++)
        {
            outmat->data[i][j] = (*inmat)->data[i][j];
        }
    }
    /* DOES NOT CLEAN */
    return outmat;
}

int VC_GEE_split(MATRIX **matptr, MATRIX **discptr, MATRIX **matarrptr[])
{ /* discriminator VC_GEE_vector assumed to be integer-valued dbls */
    int i, iVC_GEE_start, k, VC_GEE_start, end;
    if ((*discptr)->ncols != 1)
    {
        printf("GEE Error: VC_GEE_split: discriminator must be column vec.\nVC_GEE_split: ncols = %d", (*discptr)->ncols);
        return (EXIT_FAILURE);
    }

    k = 0;

    iVC_GEE_start = (int)(*discptr)->data[0][0];
    VC_GEE_start = 0;
    end = 0;
    for (i = 1; i <= (*discptr)->nrows; i++)
    {
        if (i == (*discptr)->nrows || (*discptr)->data[i][0] != iVC_GEE_start)
        {
            MATRIX *ttemp = VC_GEE_extract_rows(matptr, VC_GEE_start, end);
            (*matarrptr)[k] = VC_GEE_matcopy(&ttemp);
            VC_GEE_destroy_matrix(&ttemp);
            make_permanent((*matarrptr)[k]);
            k++;
            VC_GEE_start = end + 1;
            if (i < (*discptr)->nrows) /* don't need iVC_GEE_start at end of loop */
                iVC_GEE_start = (int)(*discptr)->data[i][0];
        }
        if (VC_GEE_start < (*discptr)->nrows)
            end++;
    }
    /* DOES NOT CLEAN */
    return k;
}

/*mask socket with VC_GEE_plugm from position (row, col)*/
int VC_GEE_plug(MATRIX **VC_GEE_plugm, MATRIX **socket, int row, int col)
{
    int pcol, prow, i, j;

    pcol = (*VC_GEE_plugm)->ncols;
    prow = (*VC_GEE_plugm)->nrows;

    if (pcol + col > (*socket)->ncols || prow + row > (*socket)->nrows)
    {
        printf("GEE Error: M+-: VC_GEE_plug: socket too small");
        return (EXIT_FAILURE);
    }

    for (i = 0; i < prow; i++)
    {
        for (j = 0; j < pcol; j++)
        {
            (*socket)->data[i + row][j + col] = (*VC_GEE_plugm)->data[i][j];
        }
    }

    free_if_ephemeral(VC_GEE_plugm);
}

MATRIX *VC_GEE_form_diag(MATRIX **vec)
{
    MATRIX *tmp;
    int i, ord;

    ord = (*vec)->nrows;
    tmp = VC_GEE_create_matrix(ord, ord, EPHEMERAL);
    for (i = 0; i < ord; i++)
        tmp->data[i][i] = (*vec)->data[i][0];
    free_if_ephemeral(vec);
    return tmp;
}

MATRIX *VC_GEE_band(MATRIX **in, int wid)
{
    MATRIX *tmp;
    int i, j;
    tmp = VC_GEE_matcopy(in);
    for (i = 0; i < (*in)->nrows; i++)
    {
        for (j = i + wid; j < (*in)->ncols; j++)
        {
            tmp->data[i][j] = (double)0.;
            if ((i < (*in)->ncols) && (j < (*in)->nrows))
            {
                tmp->data[j][i] = (double)0.;
            }
        }
    }
    free_if_ephemeral(in);
    return tmp;
}

MATRIX *VC_GEE_toeplitz(MATRIX **in)
{
    MATRIX *toep = NULL, *tin = NULL, *tmp = NULL, *ttmp = NULL;
    int n, p, inrows, incols, i, j;

    inrows = (*in)->nrows;
    incols = (*in)->ncols;

    if ((inrows > incols) ? inrows % incols : incols % inrows)
    {
        printf("GEE Error: M+-:VC_GEE_toeplitz: argument invalid");
        return (NULL);
    }

    if (inrows > incols)
    {
        p = incols;
        n = inrows / p;
        tin = VC_GEE_matcopy(in);
        free_if_ephemeral(in);
    }
    else
    {
        p = inrows;
        n = incols / p;
        tin = VC_GEE_transp(in);
    }

    toep = VC_GEE_create_matrix(n * p, n * p, EPHEMERAL);

    for (i = 0; i < n; i++)
    {
        tmp = VC_GEE_extract_rows(&tin, i * p, (i * p) + p - 1);
        make_permanent(tmp);
        ttmp = VC_GEE_transp(&tmp);
        make_permanent(ttmp);
        if (i == 0)
        {
            for (j = 0; j < n; j++)
            {
                if (inrows > incols)
                    VC_GEE_plug(&tmp, &toep, j * p, j * p);
                else
                    VC_GEE_plug(&ttmp, &toep, j * p, j * p);
            }
        }
        else
        {
            for (j = 0; j < n - i; j++)
            {
                VC_GEE_plug(&ttmp, &toep, j * p, (j + i) * p);
                VC_GEE_plug(&tmp, &toep, (j + i) * p, j * p);
            }
        }
        VC_GEE_destroy_matrix(&tmp);
        VC_GEE_destroy_matrix(&ttmp);
    }
    VC_GEE_destroy_matrix(&tin);
    return toep;
}

#define get_nelem(x) (((x)->nrows) * ((x)->ncols))

double VC_GEE_elsum(MATRIX **x)
{
    double t = 0.;
    int i, j;
    for (i = 0; i < (*x)->nrows; i++)
    {
        for (j = 0; j < (*x)->ncols; j++)
        {
            t += (*x)->data[i][j];
        }
    }
    free_if_ephemeral(x);
    return t;
}

MATRIX *VC_GEE_matabs(MATRIX **x)
{
    MATRIX *tmp;
    int i, j;
    tmp = VC_GEE_create_matrix((*x)->nrows, (*x)->ncols, EPHEMERAL);
    for (i = 0; i < (*x)->nrows; i++)
    {
        for (j = 0; j < (*x)->ncols; j++)
        {
            tmp->data[i][j] = fabs((*x)->data[i][j]);
        }
    }
    free_if_ephemeral(x);
    return tmp;
}

double VC_GEE_matmax(MATRIX **x)
{
    double t;
    int i, j;

    t = (*x)->data[0][0];
    for (i = 0; i < (*x)->nrows; i++)
    {
        for (j = 0; j < (*x)->ncols; j++)
        {
            if ((*x)->data[i][j] > t)
                t = (*x)->data[i][j];
        }
    }
    free_if_ephemeral(x);
    return t;
}

MATRIX *VC_GEE_matexp(MATRIX **x)
{
    MATRIX *tmp = VC_GEE_create_matrix((*x)->nrows, (*x)->ncols, EPHEMERAL);
    int i, j;
    for (i = 0; i < (*x)->nrows; i++)
    {
        for (j = 0; j < (*x)->ncols; j++)
        {
            tmp->data[i][j] = exp((*x)->data[i][j]);
        }
    }
    free_if_ephemeral(x);
    return tmp;
}

MATRIX *VC_GEE_matadd(MATRIX **mat1, MATRIX **mat2)
{
    int i, j;
    if (((*mat1)->ncols != (*mat2)->ncols) || ((*mat1)->nrows != (*mat2)->nrows))
    {
        printf("GEE Error: VC_GEE_matadd: args (%dx%d) + (%dx%d) don't conform.\nfatal error",
               (*mat1)->nrows, (*mat1)->ncols, (*mat2)->nrows, (*mat2)->ncols);
        return (NULL);
    }
    MATRIX *result = VC_GEE_create_matrix((*mat1)->nrows, (*mat1)->ncols, EPHEMERAL);
    for (i = 0; i < result->nrows; i++)
    {
        for (j = 0; j < result->ncols; j++)
        {
            result->data[i][j] = (*mat1)->data[i][j] + (*mat2)->data[i][j];
        }
    }
    free_if_ephemeral(mat1);
    free_if_ephemeral(mat2);
    return result;
}

MATRIX *VC_GEE_matsub(MATRIX **mat1, MATRIX **mat2)
{
    int i, j;
    if (((*mat1)->ncols != (*mat2)->ncols) || ((*mat1)->nrows != (*mat2)->nrows))
    {
        printf("GEE Error: VC_GEE_matadd: args (%dx%d) + (%dx%d) don't conform.\nfatal error",
               (*mat1)->nrows, (*mat1)->ncols, (*mat2)->nrows, (*mat2)->ncols);
        return (NULL);
    }
    MATRIX *result = VC_GEE_create_matrix((*mat1)->nrows, (*mat1)->ncols, EPHEMERAL);
    for (i = 0; i < result->nrows; i++)
    {
        for (j = 0; j < result->ncols; j++)
        {
            result->data[i][j] = (*mat1)->data[i][j] - (*mat2)->data[i][j];
        }
    }
    free_if_ephemeral(mat1);
    free_if_ephemeral(mat2);
    return result;
}

MATRIX *VC_GEE_matmult(MATRIX **mat1, MATRIX **mat2)
{
    int i, j, k;

    if ((*mat1)->ncols != (*mat2)->nrows)
    {
        printf("GEE Error: VC_GEE_matmult: args (%dx%d) * (%dx%d) don't conform.\n",
               (*mat1)->nrows, (*mat1)->ncols, (*mat2)->nrows, (*mat2)->ncols);
        return (NULL);
    }

    MATRIX *result = VC_GEE_create_matrix((*mat1)->nrows, (*mat2)->ncols, EPHEMERAL);

    for (i = 0; i < result->nrows; i++)
    {
        for (j = 0; j < result->ncols; j++)
        {
            for (k = 0; k < (*mat2)->nrows; k++)
            {
                // printf("i: %d, j: %d, k: %d\n", i, j, k);
                result->data[i][j] += (*mat1)->data[i][k] * (*mat2)->data[k][j];
            }
        }
    }
    free_if_ephemeral(mat1);
    free_if_ephemeral(mat2);
    return result;
}

MATRIX *VC_GEE_px1_times_pxq(MATRIX **px1, MATRIX **pxq) /* mult elements of a colvec */
/* across corresp row of mat */
{
    MATRIX *tmp;
    double *load, colel;
    int i, j;

    if ((*px1)->ncols != 1)
    {
        printf("GEE Error: M+-: VC_GEE_px1_times_pxq: arg1 not a col-vec");
        return (NULL);
    }
    if ((*px1)->nrows != (*pxq)->nrows)
    {
        printf("GEE Error: M+-: VC_GEE_px1_times_pxq: args not conforming");
        return (NULL);
    }
    tmp = VC_GEE_matcopy((pxq));
    for (i = 0; i < tmp->nrows; i++)
    {
        colel = (*px1)->data[i][0];
        for (j = 0; j < tmp->ncols; j++)
        {
            tmp->data[i][j] *= colel;
        }
    }
    free_if_ephemeral(px1);
    free_if_ephemeral(pxq);
    return tmp;
}

MATRIX *VC_GEE_pxq_divby_px1(MATRIX **pxq, MATRIX **px1) /* divide elements of a colvec */
/* into corresp row of mat */
{
    MATRIX *tmp;
    double *load, colel;
    int i, j;
    if ((*px1)->ncols != 1)
    {
        printf("GEE Error: M+-: VC_GEE_pxq_divby_px1: arg2 not a col-vec");
        return (NULL);
    }
    if ((*px1)->nrows != (*pxq)->nrows)
    {
        printf("GEE Error: M+-: VC_GEE_pxq_divby_px1: args not conforming");
        return (NULL);
    }

    tmp = VC_GEE_matcopy(pxq);
    for (i = 0; i < tmp->nrows; i++)
    {
        colel = (*px1)->data[i][0];
        for (j = 0; j < tmp->ncols; j++)
        {
            tmp->data[i][j] /= colel;
        }
    }
    free_if_ephemeral(px1);
    free_if_ephemeral(pxq);
    return tmp;
}

MATRIX *VC_GEE_scalar_times_matrix(double a, MATRIX **X)
{

    int i, j;
    MATRIX *tmp = VC_GEE_matcopy(X);
    for (i = 0; i < tmp->nrows; i++)
    {
        for (j = 0; j < tmp->ncols; j++)
        {
            tmp->data[i][j] *= a;
        }
    }
    free_if_ephemeral(X);
    return tmp;
}

int VC_GEE_matdump(MATRIX **mat)
{
    if ((*mat) == NULL)
    {
        printf("warning: VC_GEE_matdump: arg is NULL");
        return (EXIT_FAILURE);
    }

    for (int i = 0; i < (*mat)->nrows; i++)
    {
        for (int j = 0; j < (*mat)->ncols; j++)
        {
            printf("%.3f ", (*mat)->data[i][j]);
        }
        printf("\n");
    }
    return (EXIT_SUCCESS);
    /* DOES NOT CLEAN */
}

MATRIX *VC_GEE_covlag(MATRIX **inmat, int lag, int demean)
{
    MATRIX *xrows[MAX_COVLAG], *res, *temp, *temp1, *temp2;
    int n, i, j, nv, q;
    double nrec;

    n = (*inmat)->nrows;
    nrec = (double)1. / (double)n;
    if (n > MAX_COVLAG)
    {
        printf("GEE Error: VC_GEE_covlag: arg has > MAX_COVLAG rows");
        return (NULL);
    }

    nv = (*inmat)->ncols;

    res = VC_GEE_create_matrix(nv, lag * nv, EPHEMERAL);

    for (q = 0; q < n; q++)
    {
        xrows[q] = VC_GEE_extract_rows(inmat, q, q);
        make_permanent(xrows[q]);
    }

    for (i = 0; i < lag; i++)
    {
        temp = VC_GEE_create_matrix(nv, nv, EPHEMERAL);
        for (j = i; j < n; j++)
        {
            if ((j - i) < n)
            {
                temp1 = VC_GEE_transp(&xrows[j]);
                temp2 = VC_GEE_matmult(&temp1, &xrows[j - i]);
                temp = VC_GEE_matadd(&temp, &temp2);

                // temp = VC_GEE_matadd(
                //         temp,
                //         VC_GEE_matmult(VC_GEE_transp(xrows[j]), xrows[j - i]));
            }
        }

        temp1 = VC_GEE_scalar_times_matrix(nrec, &temp);
        VC_GEE_plug(&temp1, &res, 0, i * nv);
    }

    for (q = 0; q < n; q++)
    {
        VC_GEE_destroy_matrix(&xrows[q]);
    }
    return res;
}

MATRIX *VC_GEE_ident(int ord)
{
    MATRIX *I;
    int i;

    I = VC_GEE_create_matrix(ord, ord, EPHEMERAL);
    for (i = 0; i < ord; i++)
        I->data[i][i] = (double)1.0;
    return I;
}

MATRIX *VC_GEE_col_1s(int k)
{
    MATRIX *tmp;
    int i;
    tmp = VC_GEE_create_matrix(k, 1, EPHEMERAL);
    for (i = 0; i < k; i++)
    {
        tmp->data[i][0] = 1.;
    }
    return tmp;
}

int VC_GEE_nchanges(MATRIX **X)
{
    /* returns integer telling how often the value of X */
    /* changes from row to row.  X must be column VC_GEE_vector */

    int tmp = 1, iVC_GEE_start, i;

    if ((*X)->ncols != 1)
    {
        printf("GEE Error: VC_GEE_nchanges:  must be column VC_GEE_vector; ncols = %d",
               (*X)->ncols);
        return (EXIT_FAILURE);
    }

    iVC_GEE_start = (int)(*X)->data[0][0];

    for (i = 1; i < (*X)->nrows; i++)
    {
        if ((*X)->data[i][0] != iVC_GEE_start)
        {
            tmp++;
            iVC_GEE_start = (int)(*X)->data[i][0];
        }
    }
    return tmp;
}

MATRIX *VC_GEE_matanticlog(MATRIX **x)
{
    int i, j;
    MATRIX *tmp = VC_GEE_create_matrix((*x)->nrows, (*x)->ncols, EPHEMERAL);
    for (i = 0; i < (*x)->nrows; i++)
    {
        for (j = 0; j < (*x)->ncols; j++)
        {
            tmp->data[i][j] = 1 - exp(-exp((*x)->data[i][j]));
        }
    }
    free_if_ephemeral(x);
    return tmp;
}

/* gee support @(#) VC_GEE_diag_as_vec.c 3.3 94/03/09 */

MATRIX *VC_GEE_diag_as_vec(MATRIX **inmat)
{
    int i;
    MATRIX *outmat;

    if ((*inmat)->ncols != (*inmat)->nrows)
    {
        printf("GEE Error: M+-: VC_GEE_diag_as_vec: arg is not a square matrix");
        return (NULL);
    }

    outmat = VC_GEE_create_matrix((*inmat)->nrows, 1, EPHEMERAL);
    for (i = 0; i < (*inmat)->nrows; i++)
    {
        outmat->data[i][0] = (*inmat)->data[i][i];
    }
    free_if_ephemeral(inmat);
    return outmat;
}

MATRIX *VC_GEE_matsqrt(MATRIX **x)
{
    int i, j;
    MATRIX *tmp;
    tmp = VC_GEE_matcopy(x);
    for (i = 0; i < (*x)->nrows; i++)
    {
        for (j = 0; j < (*x)->ncols; j++)
        {
            tmp->data[i][j] = sqrt((*x)->data[i][j]);
        }
    }
    free_if_ephemeral(x);
    return tmp;
}

MATRIX *VC_GEE_mat1over(MATRIX **x)
{
    int i, j;
    MATRIX *tmp = VC_GEE_matcopy(x);
    for (i = 0; i < (*x)->nrows; i++)
    {
        for (j = 0; j < (*x)->ncols; j++)
        {
            tmp->data[i][j] = 1. / ((*x)->data[i][j]);
        }
    }
    free_if_ephemeral(x);
    return tmp;
}

MATRIX *VC_GEE_matnpdf(MATRIX **x)
{
    int i, j;
    MATRIX *tmp = VC_GEE_create_matrix((*x)->nrows, (*x)->ncols, EPHEMERAL);
    for (i = 0; i < (*x)->nrows; i++)
    {
        for (j = 0; j < (*x)->ncols; j++)
        {
            tmp->data[i][j] = gsl_ran_ugaussian_pdf((*x)->data[i][j]);
        }
    }
    free_if_ephemeral(x);
    return tmp;
}

MATRIX *VC_GEE_matncdf(MATRIX **x)
{
    int i, j;
    MATRIX *tmp = VC_GEE_create_matrix((*x)->nrows, (*x)->ncols, EPHEMERAL);
    for (i = 0; i < (*x)->nrows; i++)
    {
        for (j = 0; j < (*x)->ncols; j++)
        {
            tmp->data[i][j] = gsl_cdf_ugaussian_P((*x)->data[i][j]);
        }
    }
    free_if_ephemeral(x);
    return tmp;
}

void VC_GEE_destroy_matrix(MATRIX **mat)
{
    if (*mat == NULL)
    {
        return;
    }
    for (int i = 0; i < (*mat)->nrows; i++)
    {
        Free((*mat)->data[i]);
    }
    Free((*mat)->data);

    Free((*mat));
    *mat = NULL;
}

void VC_GEE_destroy_double_matrix(MATRIX ***mat, int nmat)
{
    for (int i = 0; i < nmat; i++)
    {
        VC_GEE_destroy_matrix(&((*mat)[i]));
    }
    Free(*mat);
    *mat = NULL;
}

void VC_GEE_destroy_matrices(int nmat, ...)
{
    va_list ap;
    va_start(ap, nmat);
    for (int i = 0; i < nmat; i++)
    {
        MATRIX **mat = va_arg(ap, MATRIX **);
        VC_GEE_destroy_matrix(mat);
    }
    va_end(ap);
}

void free_if_ephemeral(MATRIX **mat)
{
    if (is_ephemeral(*mat))
    {
        VC_GEE_destroy_matrix(mat);
    }
}

// Custom error handler function
void custom_error_handler(const char *reason, const char *file, int line, int gsl_errno)
{
    if (_GSL_ERROR_FLAG)
    {
        printf("GSL Error: %s\n", reason);
        fprintf(stderr, "Location: %s:%d\n", file, line);
        fprintf(stderr, "GSL Error Code: %d\n", gsl_errno);
    }
    // You can add more custom error handling here if needed
}

// Function to convert MATRIX to gsl_matrix
gsl_matrix *matrix_to_gsl(MATRIX **mat)
{
    gsl_matrix *gslMat = gsl_matrix_alloc((*mat)->nrows, (*mat)->ncols);
    for (int i = 0; i < (*mat)->nrows; i++)
    {
        for (int j = 0; j < (*mat)->ncols; j++)
        {
            gsl_matrix_set(gslMat, i, j, (*mat)->data[i][j]);
        }
    }
    return gslMat;
}

MATRIX *gsl_to_matrix(gsl_matrix *gslMat)
{
    MATRIX *mat = VC_GEE_create_matrix(gslMat->size1, gslMat->size2, EPHEMERAL);
    for (int i = 0; i < mat->nrows; i++)
    {
        for (int j = 0; j < mat->ncols; j++)
        {
            mat->data[i][j] = gsl_matrix_get(gslMat, i, j);
        }
    }
    return mat;
}

MATRIX *matrix_inverse(MATRIX **mat)
{
    // Set the custom error handler so that it won't abort the program
    gsl_set_error_handler(&custom_error_handler);

    gsl_matrix *gslMat = matrix_to_gsl(mat);

    // Perform LU decomposition
    gsl_permutation *p = gsl_permutation_alloc((*mat)->nrows);
    int signum;
    gsl_linalg_LU_decomp(gslMat, p, &signum);

    // Invert the matrix
    gsl_matrix *inverseGSL = gsl_matrix_alloc((*mat)->nrows, (*mat)->ncols);
    int status = gsl_linalg_LU_invert(gslMat, p, inverseGSL);
    // Check for inversion failure
    if (status != GSL_SUCCESS)
    {
        if (_GSL_ERROR_FLAG)
            printf("Error: Matrix inversion failed (Status: %d)\n", status);
        gsl_matrix_free(gslMat);
        gsl_matrix_free(inverseGSL);
        gsl_permutation_free(p);
        return NULL; // or handle the error appropriately
    }

    // Convert GSL matrix back to MATRIX
    MATRIX *inverseMat = gsl_to_matrix(inverseGSL);

    // Free GSL resources
    gsl_matrix_free(gslMat);
    gsl_matrix_free(inverseGSL);
    gsl_permutation_free(p);

    free_if_ephemeral(mat);

    return inverseMat;
}

void PrintMatrix(MATRIX **mat)
{
    if ((*mat) == NULL)
    {
        printf("Matrix is NULL\n");
        return;
    }
    printf("Matrix(%d x %d):\n", (*mat)->nrows, (*mat)->ncols);
    for (int i = 0; i < (*mat)->nrows && i < 20; i++)
    {
        for (int j = 0; j < (*mat)->ncols && j < 20; j++)
        {
            printf("%.3f ", (*mat)->data[i][j]);
        }
        printf("\n");
    }
}
