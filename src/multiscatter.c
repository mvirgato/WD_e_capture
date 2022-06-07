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
#include "multiscatter.h"
#include "cap_funcs.h"
#include "crate_apprx_funcs.h"
#include "cap_interp.h"

//================================================================
//FUNCTIONS
//================================================================

double nstar(double B, double muf, double mchi, int oper)
{
    return 1./(1. - exp(-mstarInterp(B, muf, oper)/mchi ) );
}

//================================================================
// INTEGRANDS
//================================================================

// double simple_integrand(){

// }

double monteIntegrandMultiScatt(double x[], size_t dim, void *p)
{
    (void)(dim);

    struct crate_params *params = (struct crate_params *)p;

    double mdm = (params->mchi);
    int np = (params->npts);
    int oper = (params->oper);
    double z0 = (params->z0);
    double k0 = (params->k0);

    double r = x[0]; // r/Rstar
    // double u = x[1]; // uchi
    double E = x[1]; // Ep = (Eumax - Eumin)*Eu + Eumin
    double s = x[2]; // sp = (smax - smin)*s + smin
    double t = x[3]; // tp = -tmin*t + tmin

    double muf = muFeinterp(r, np);
    double ne = nEinterp(r, np);
    double B = Binterp(r, np);

    double monteInt = 0.5*integrandUappxNoPB(r, E, s, t, muf, ne, B, mdm, oper)/nstar(B, muf, mdm, oper);

    return monteInt;
}

//================================================================
// INTEGRATORS
//================================================================

double crateMultiScatt(double mchi, int oper, int npts)
{

    double res, err;

    struct crate_params params = {mchi, oper, npts};

    // double term_max[5] = {1, umax(), 1, 1, 1};
    // double term_min[5] = {0, 0, 0, 0, 0};

    // const gsl_rng_type *Trng;
    // gsl_rng *rng;

    // gsl_monte_function multiFunc = {&monteIntegrandMultiScatt, 5, &params};

    // size_t calls = 10000;

    // Trng = gsl_rng_default;
    // rng = gsl_rng_alloc(Trng);

    // {
    //     gsl_monte_plain_state *s = gsl_monte_plain_alloc(5);
    //     gsl_monte_plain_integrate(&multiFunc, term_min, term_max, 5, calls, rng, s, &res, &err);
    //     gsl_monte_plain_free(s);

    //     // display_results("plain", res, err);
    // }

    double term_max[4] = {1, 1, 1, 1};
    double term_min[4] = { 0, 0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function multiFunc = {&monteIntegrandMultiScatt, 4, &params};

    size_t calls = 10000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    {
        gsl_monte_plain_state *s = gsl_monte_plain_alloc(4);
        gsl_monte_plain_integrate(&multiFunc, term_min, term_max, 4, calls, rng, s, &res, &err);
        gsl_monte_plain_free(s);

        // display_results("plain", res, err);
    }

    // {
    //     gsl_monte_miser_state *s = gsl_monte_miser_alloc(5);
    //     gsl_monte_miser_integrate(&multiFunc, term_min, term_max, 5, calls, rng, s, &res, &err);
    //     gsl_monte_miser_free(s);

    //     display_results("misner", res, err);
    // }

    // {
    //     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(4);

    //     gsl_monte_vegas_integrate(&multiFunc, term_min, term_max, 4, calls / 5, rng, s, &res, &err);
    //     display_results("vegas warm-up", res, err);
    //     printf("converging...\n");

    //     do
    //     {
    //         gsl_monte_vegas_integrate(&multiFunc, term_min, term_max, 4, calls, rng, s, &res, &err);
    //         printf("result = % .6e sigma = % .6e "
    //                "chisq/dof = %.1e\n",
    //                res, err, gsl_monte_vegas_chisq(s));
    //     }

    //     while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.2);

    //     display_results("vegas final", res, err);

    //     gsl_monte_vegas_free(s);
    // }

    gsl_rng_free(rng);

    return res* CoeffUappx(mchi);
}