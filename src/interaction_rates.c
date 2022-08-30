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
#include "interaction_rates.h"
#include "collective_effects.h"

double intRatePrefac(double mchi){

    return 0.5 * mE * mE / mchi / pi /pi ;
}

double integrandIntRate(double r, double Ep, double sp, double tp, double muFe, double nE, double B, double mchi, int oper, double k0, double z0)
{
    // Full 3D integration

    double E1 = Eumin(muFe, mchi, k0);
    double E2 = Eumax(muFe);
    double Eu = (E2 - E1) * Ep + E1;

    double s1 = smin(Eu, 0, B, mchi); // smin
    double s2 = smax(Eu, 0, B, mchi); // smax
    double s  = (s2 - s1) * sp + s1;

    double t1 = (1. - fmax(1. - z0 * mchi / mt, 0))*tmin(s, mchi);
    double t2 = tmax();
    double t  = -t1 * tp + t1;

    double EstJac = (E2 - E1) * (s2 - s1) * (t2 - t1); // Jacobian to change (Eu, s, t) to unit cube
    // printf("%0.5e\n",tp);
    double soln = solnFunc(Eu, s, t, 0, mchi, B, muFe);
    // * HeavisidePhaseSpace(Eu, s, B, s1, s2)
           
    double int_main = CollEffectsFF(t, 0, mchi, B, muFe, nE, soln) * (zeta(nE, muFe) * Eu * s * opersdXS(s, t, mchi, oper) * HeavisidePB(Eu, s, t, 0, mchi, B, muFe, soln) * sqrt(B / ( 1.0 - B ) )/ beta_cap(s, mchi) / gamma_cap(s, mchi));

    // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\n", Eu, B, s1);

    // printf("%f\n", r);
    return EstJac * int_main;
}

double monteIntegrandIntRate(double x[], size_t dim, void *p)
{
    (void)(dim);

    struct int_rate_params *params = (struct int_rate_params *)p;

    double r   = (params->r); // r/Rstar
    double mdm = (params->mchi);
    int np     = (params->npts);
    int oper   = (params->oper);
    double z0  = (params->z0);
    double k0  = (params->k0);

    double E = x[0]; // Ep = (Eumax - Eumin)*Eu + Eumin
    double s = x[1]; // sp = (smax - smin)*s + smin
    double t = x[2]; // tp = -tmin*t + tmin

    double muf = muFeinterp(r, np);
    double ne  = nEinterp(r, np);
    double B   = Binterp(r, np);

    double monteInt = integrandIntRate(r, E, s, t, muf, ne, B, mdm, oper, k0, z0);

    return monteInt;
}

double intRate(double r, double mchi, int oper, int npts, void *cont_vars)
{
    // Input r in r/Rstar

    double res, err;

    struct cont_pars *contp  = (struct cont_pars *)cont_vars;

    double z0 = (contp->z0);
    double k0 = (contp->k0);

    struct int_rate_params params = {r, mchi, oper, npts, z0, k0};

    double term_max[3] = {1, 1, 1};
    double term_min[3] = {0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function IntRateIntegrand = {&monteIntegrandIntRate, 3, &params};

    size_t calls = 50000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    // {
    //     gsl_monte_plain_state *s = gsl_monte_plain_alloc(3);
    //     gsl_monte_plain_integrate(&IntRateIntegrand, term_min, term_max, 3, calls, rng, s, &res, &err);
    //     gsl_monte_plain_free(s);

    //     // display_results("plain", res, err);
    // }

    // {
    //     gsl_monte_miser_state *s = gsl_monte_miser_alloc(3);
    //     gsl_monte_miser_integrate(&IntRateIntegrand, term_min, term_max, 3, calls, rng, s, &res, &err);
    //     gsl_monte_miser_free(s);

    //     // display_results("misner", res, err);
    //     // printf("MISER RESULT:\t%0.5e +- %0.3e\n", res*CoeffUappx(mchi), err*CoeffUappx(mchi));
    // }

    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(3);

        gsl_monte_vegas_integrate(&IntRateIntegrand, term_min, term_max, 3, calls / 5, rng, s, &res, &err);
        // display_results("vegas warm-up", res, err);
        // printf("converging...\n");

        do
        {
            gsl_monte_vegas_integrate(&IntRateIntegrand, term_min, term_max, 3, calls, rng, s, &res, &err);
            // printf("result = % .6e sigma = % .6e "
            //        "chisq/dof = %.1e\n",
            //        res, err, gsl_monte_vegas_chisq(s));
        }

        while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);

        // display_results("vegas final", res, err);

        gsl_monte_vegas_free(s);
        printf("VEGAS RESULT:\t%0.5e +- %0.3e\n", res, err);
    }

    gsl_rng_free(rng);

    return res * intRatePrefac(mchi);
}


