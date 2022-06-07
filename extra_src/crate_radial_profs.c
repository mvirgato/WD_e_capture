#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>

#include "globals.h"
#include "cap_funcs.h"
#include "crate_apprx_funcs.h"
#include "cap_interp.h"

//=====================================================
// MONTE RADIAL PROFILE INTEGRANDS
//=====================================================

double monteIntegrandCrateRadFull(double x[], size_t dim, void *p)
{
    (void)(dim);

    struct radprof_params *params = (struct radprof_params *)p;

    double r   = (params->rad);
    double mdm = (params->mchi);
    int np     = (params->npts);
    int oper   = (params->oper);

    double E = x[0]; // Ep = (Eumax - Eumin)*Eu + Eumin
    double s = x[1]; // sp = (smax - smin)*s + smin
    double t = x[2]; // tp = -tmin*t + tmin

    double muf = muFeinterp(r, np);
    double ne  = nEinterp(r, np);
    double B   = Binterp(r, np);

    double monteInt = integrandUappx(r, E, s, t, muf, ne, B, mdm, oper);

    return monteInt;
}

double monteIntegrandCrateRadNoPB(double x[], size_t dim, void *p)
{
    (void)(dim);

    struct radprof_params *params = (struct radprof_params *)p;

    double r = (params->rad);
    double mdm = (params->mchi);
    int np = (params->npts);
    int oper = (params->oper);

    double E = x[0]; // Ep = (Eumax - Eumin)*Eu + Eumin
    double s = x[1]; // sp = (smax - smin)*s + smin
    double t = x[2]; // tp = -tmin*t + tmin

    double muf = muFeinterp(r, np);
    double ne = nEinterp(r, np);
    double B = Binterp(r, np);

    double monteInt = integrandNoPB(r, 0, E, s, t, muf, ne, B, mdm, oper);

    return monteInt;
}

double monteIntegrandApprx(double x[], size_t dim, void *p)
{
    (void)(dim);

    struct radprof_params *params = (struct radprof_params *)p;

    double r = (params->rad);
    double mdm = (params->mchi);
    int np = (params->npts);
    int oper = (params->oper);

    double u = x[0]; // uchi
    double E = x[1]; // Ep = (Eumax - Eumin)*Eu + Eumin
    double t = x[2]; // tp = -tmin*t + tmin

    double muf = muFeinterp(r, np);
    double ne = nEinterp(r, np);
    double B = Binterp(r, np);

    double monteInt = integrandApprx(r, u, E, t, muf, ne, B, mdm, oper);

    return monteInt;
}

//////////////////////////////////////////////////////////
// RADIAL PROFILE INTEGRATORS
//////////////////////////////////////////////////////////

double crateRaProfFull(double rad, double mchi, int oper, int npts)
{

    double res, err;

    struct radprof_params params = {rad, mchi, oper, npts};

    double term_max[3] = {1, 1, 1};
    double term_min[3] = {0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function cradFullMfunc = {&monteIntegrandCrateRadFull, 3, &params};

    size_t calls = 50000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    // {
    //     gsl_monte_plain_state *s = gsl_monte_plain_alloc(4);
    //     gsl_monte_plain_integrate(&cradFullMfunc, term_min, term_max, 4, calls, rng, s, &res, &err);
    //     gsl_monte_plain_free(s);

    //     display_results("plain", res, err);
    // }

    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(3);
        gsl_monte_miser_integrate(&cradFullMfunc, term_min, term_max, 3, calls, rng, s, &res, &err);
        gsl_monte_miser_free(s);

        // display_results("misner", res, err);
    }

    // {
    //     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(4);

    //     gsl_monte_vegas_integrate(&cradFullMfunc, term_min, term_max, 4, calls / 5, rng, s, &res, &err);
    //     // display_results("vegas warm-up", res, err);
    //     printf("converging...\n");
    //     printf("initial iters: %d\n", iters);

    //     int MAX_ITERS = 20;

    //     do
    //     {
    //         gsl_monte_vegas_integrate(&cradFullMfunc, term_min, term_max, 4, calls, rng, s, &res, &err);
    //         // printf("result = % .6e sigma = % .6e "
    //         //        "chisq/dof = %.1e\n",
    //         //        res, err, gsl_monte_vegas_chisq(s));

    //         iters++;

    //     }

    //     while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.2 && iters < MAX_ITERS);

    //     if (iters >= MAX_ITERS)
    //     {
    //         printf("Max iterations reached without convergence\n");
    //     }
    //     printf("final iters: %d\n", iters);

    //     display_results("vegas final", res, err);

    //     gsl_monte_vegas_free(s);
    // }

    gsl_rng_free(rng);

    return res * Coeff(mchi);
}

double crateRaProfNoPb(double rad, double mchi, int oper, int npts)
{

    double res, err;

    struct radprof_params params = {rad, mchi, oper, npts};

    double term_max[4] = {umax(), 1, 1, 1};
    double term_min[4] = {0, 0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function cradFullMfunc = {&monteIntegrandCrateRadNoPB, 4, &params};

    size_t calls = 10000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    {
        gsl_monte_plain_state *s = gsl_monte_plain_alloc(4);
        gsl_monte_plain_integrate(&cradFullMfunc, term_min, term_max, 4, calls, rng, s, &res, &err);
        gsl_monte_plain_free(s);

        display_results("plain", res, err);
    }

    // {
    //     gsl_monte_miser_state *s = gsl_monte_miser_alloc(4);
    //     gsl_monte_miser_integrate(&cradFullMfunc, term_min, term_max, 4, calls, rng, s, &res, &err);
    //     gsl_monte_miser_free(s);

    //     display_results("misner", res, err);
    // }

    // {
    //     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(4);

    //     gsl_monte_vegas_integrate(&cradFullMfunc, term_min, term_max, 4, calls / 5, rng, s, &res, &err);
    //     // display_results("vegas warm-up", res, err);
    //     printf("converging...\n");

    //     int MAX_ITERS = 50;

    //     do
    //     {
    //         gsl_monte_vegas_integrate(&cradFullMfunc, term_min, term_max, 4, calls, rng, s, &res, &err);
    //         // printf("result = % .6e sigma = % .6e "
    //         //        "chisq/dof = %.1e\n",
    //         //        res, err, gsl_monte_vegas_chisq(s));

    //         iters++;

    //     }

    //     while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.1 && iters < MAX_ITERS);

    //     if (iters >= MAX_ITERS){
    //         printf("Max iterations reached without convergence\n");
    //     }
    //     printf("%d", iters);

    //     display_results("vegas final", res, err);

    //     gsl_monte_vegas_free(s);
    // }

    gsl_rng_free(rng);

    return res * Coeff(mchi);
}

double crateRaProfAppx(double rad, double mchi, int oper, int npts, int iters)
{

    struct radprof_params params = {rad, mchi, oper, npts};

    double res, err;

    double term_max[3] = {umax(), 1, 1};
    double term_min[3] = {0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function appxfunc = {&monteIntegrandApprx, 3, &params};

    size_t calls = 10000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    // {
    //     gsl_monte_plain_state *s = gsl_monte_plain_alloc(3);
    //     gsl_monte_plain_integrate(&appxfunc, term_min, term_max, 3, calls, rng, s, &res, &err);
    //     gsl_monte_plain_free(s);

    //     display_results("plain", res, err);
    // }

    // {
    //     gsl_monte_miser_state *s = gsl_monte_miser_alloc(3);
    //     gsl_monte_miser_integrate(&appxfunc, term_min, term_max, 3, calls, rng, s, &res, &err);
    //     gsl_monte_miser_free(s);

    //     display_results("misner", res, err);
    // }

    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(3);

        gsl_monte_vegas_integrate(&appxfunc, term_min, term_max, 3, calls / 5, rng, s, &res, &err);
        // display_results("vegas warm-up", res, err);
        printf("converging...\n");

        int MAX_ITERS = 50;

        do
        {
            gsl_monte_vegas_integrate(&appxfunc, term_min, term_max, 3, calls, rng, s, &res, &err);
            // printf("result = % .6e sigma = % .6e "
            //        "chisq/dof = %.1e\n",
            //        res, err, gsl_monte_vegas_chisq(s));

            iters++;

        }

        while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.2 && iters < MAX_ITERS);

        if (iters >= MAX_ITERS)
        {
            printf("Max iterations reached without convergence\n");
        }
        printf("final iters: %d\n", iters);

        display_results("vegas final", res, err);

        gsl_monte_vegas_free(s);
    }

    gsl_rng_free(rng);

    return res * Coeff(mchi);
}
