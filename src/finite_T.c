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
#include "finite_T.h"
#include "cap_funcs.h"
#include "crate_apprx_funcs.h"
#include "cap_interp.h"

//=============================================================
// Misc. Functions
//=============================================================

double Coeff_temp(double mchi, double temp)
{
    // numerical factors

    return (2. * rhoDM * sqrt(3. / 2. / pi / (vd * vd * kmTOm * kmTOm / cspeed / cspeed + 3. * temp / mt)) * (RS * RS * RS * mTOinveV * mTOinveV * mTOinveV) * mt * mt) / ((vs * kmTOm / cspeed) * (cmTOm * cmTOm * cmTOm * mTOinveV * mTOinveV * mTOinveV) * pi * mchi * mchi * inveVTOs);
}

double CoeffUappx_T(double mchi, double temp)
{
    double vs1 = vs*kmTOm/cspeed;
    double vd1 = vd*kmTOm/cspeed;
    return (2. * rhoDM * erf(sqrt(3. / 2./(vd1*vd1 +3.*temp/mt)) * vs1) * mt * mt * (RS * RS * RS * mTOinveV * mTOinveV * mTOinveV)) / (pi * (vs1) * mchi * mchi * cmTOm * cmTOm * cmTOm * mTOinveV * mTOinveV * mTOinveV * inveVTOs);
}

double fMB_temp(double uchi, double temp)
{
    // Maxwell Boltzmann distribution  x uchi without constants

    return (exp(-3. * (uchi - vs * kmTOm / cspeed) * (uchi - vs * kmTOm / cspeed) / (2. * (3. * temp / mt + vd * vd * kmTOm * kmTOm / cspeed / cspeed))) - exp(-3. * (uchi + vs * kmTOm / cspeed) * (uchi + vs * kmTOm / cspeed) / (2. * (3. * temp / mt + vd * vd * kmTOm * kmTOm / cspeed / cspeed))));
}

//=============================================================
// Integration Limits
//=============================================================

double EumaxT(double muFe, double temp, double ct){
    return 1 + muFe/mt + ct*temp/mt;
}

double EuminT(double muFe, double mchi, double temp, double k0, double ct){
     return fmax(1., 1. + (EumaxT(muFe, temp, ct) - ct * temp / mt - 1.)*(1. - k0 * mchi / mt - ct * temp / mt));
}

double umax_temp(double temp){
    return 10.*sqrt(vd*vd*kmTOm*kmTOm/cspeed/cspeed + 3.*temp/mt);
}

//=============================================================
// Fermi-Dirac distributions
//=============================================================

double FD(double energy, double temp){
    return 1./( 1. + exp( energy/temp ) );
}

double FDfacInitial(double Eu, double muFe, double temp){

    return 1./(1. + exp((Eu*mt - muFe - mt) / temp));
}

double FDfacFinal(double Eu, double s, double t, double uchi, double mchi, double B, double muFe, double temp){

    double E1 = mchi * (1 + uchi * uchi / 2.) / sqrt(B);
    double E2 = mt * Eu;

    double mchi2 = mchi * mchi;
    double mchi4 = mchi2 * mchi2;
    double mchi6 = mchi2 * mchi2 * mchi2;

    double mt2 = mt * mt;
    double mt4 = mt2 * mt2;
    double mt6 = mt2 * mt2 * mt2;

    double E12 = E1 * E1;
    double E22 = E2 * E2;

    double soln;

    if (s > mt2 + mchi2 + 2. * E1 * E2)
    {
        soln = (-(E2 * mchi4 * t) + E2 * mt4 * t - 2 * E2 * mt2 * s * t + E2 * s * s * t + 2 * E12 * E2 * (mt4 + (mchi2 - s) * (mchi2 - s - t) + mt2 * (-2 * mchi2 - 2 * s + t)) + sqrt(pow(2 * E1 * E2 + mchi2 + mt2 - s, 2) * (4 * E22 * mchi2 + mchi4 + 4 * E12 * mt2 + mt4 + 4 * E1 * E2 * (mchi2 + mt2 - s) - 2 * mchi2 * s + s * s - 2 * mt2 * (mchi2 + s)) * t * (mchi4 + mt4 - 2 * mchi2 * s - 2 * mt2 * (mchi2 + s) + s * (s + t))) + E1 * (pow(mchi, 6) + pow(mt, 6) + mt4 * (-mchi2 - 3 * s + t) - mchi4 * (3 * s + t) - mt2 * (mchi4 + 2 * mchi2 * s - 3 * s * s - 2 * E22 * t) - s * (s * s + 2 * E22 * t + s * t) + mchi2 * (3 * s * s - 2 * E22 * t + 2 * s * t))) / ((2 * E1 * E2 + mchi2 + mt2 - s) * (mt4 + pow(mchi2 - s, 2) - 2 * mt2 * (mchi2 + s)));
    }
    else
    {
        soln = (-(E2 * mchi4 * t) + E2 * mt4 * t - 2 * E2 * mt2 * s * t + E2 * s * s * t + 2 * E12 * E2 * (mt4 + (mchi2 - s) * (mchi2 - s - t) + mt2 * (-2 * mchi2 - 2 * s + t)) - sqrt(pow(2 * E1 * E2 + mchi2 + mt2 - s, 2) * (4 * E22 * mchi2 + mchi4 + 4 * E12 * mt2 + mt4 + 4 * E1 * E2 * (mchi2 + mt2 - s) - 2 * mchi2 * s + s * s - 2 * mt2 * (mchi2 + s)) * t * (mchi4 + mt4 - 2 * mchi2 * s - 2 * mt2 * (mchi2 + s) + s * (s + t))) + E1 * (pow(mchi, 6) + pow(mt, 6) + mt4 * (-mchi2 - 3 * s + t) - mchi4 * (3 * s + t) - mt2 * (mchi4 + 2 * mchi2 * s - 3 * s * s - 2 * E22 * t) - s * (s * s + 2 * E22 * t + s * t) + mchi2 * (3 * s * s - 2 * E22 * t + 2 * s * t))) / ((2 * E1 * E2 + mchi2 + mt2 - s) * (mt4 + pow(mchi2 - s, 2) - 2 * mt2 * (mchi2 + s)));
    }

    double EF = E2 + E1 - soln - mt - muFe;
    double finite_q0 = E1 - soln - 0.5 * mchi * uchi * uchi / sqrt(B);

    if (isnan(EF) == 0) 
    {
        if (finite_q0>0){
            return 1. / (1. + exp(-EF / temp));
        }
        else
        { 
            return 0;
        }
    }
    else
    {
        printf("Issue in finite T FD function at with:\nmdm = %0.5e, B = %0.5e, muFe = %0.5e\n", mchi, B, muFe);
        return 0.;
    }
}

//=====================================================
// CAPTURE INTEGRANDS
//=====================================================

double integrandNoPB_T(double r, double uchi, double Ep, double sp, double tp, double muFe, double nE, double B, double mchi, double temp, int oper, double ct)
{

    // No Pauli Blocking Factor

    double E1 = 1.;
    double E2 = EumaxT(muFe, temp, ct);
    double Eu = (E2 - E1) * Ep + E1;

    double s1 = smin(Eu, uchi, B, mchi); // smin
    double s2 = smax(Eu, uchi, B, mchi); // smax
    double s = (s2 - s1) * sp + s1;

    double t1 = tmin(s, mchi);
    double t2 = tmax();
    double t = -t1 * tp + t1;

    double EstJac = (E2 - E1) * (s2 - s1) * (t2 - t1); // Jacobian to change (s,t) to unit cube
    // printf("%0.5e\n",tp);
    double soln = solnFunc(Eu, s, t, uchi, mchi, B, muFe);

    double int_main = (r * r * zeta(nE, muFe) * fMB(uchi) * Eu * s * opersdXS(s, t, mchi, oper) * FDfacInitial( Eu, muFe, temp)* HeavisideUchi(Eu, s, t, uchi, mchi, B, muFe, soln) / beta_cap(s, mchi) / gamma_cap(s, mchi) / sqrt(B));

    // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\n", Eu, uchi, B, s1);

    return EstJac * int_main;
}

double integrand_T(double r, double uchi, double Ep, double sp, double tp, double muFe, double nE, double B, double mchi, double temp, int oper, double z0, double k0, double ct)
{
    // Withought fMB(u,T), placed in Monte Integrand

    double E1 = EuminT(muFe, mchi, temp, k0, ct);
    double E2 = EumaxT(muFe, temp, ct);
    double Eu = (E2 - E1) * Ep + E1;

    double s1 = smin(Eu, uchi, B, mchi); // smin
    double s2 = smax(Eu, uchi, B, mchi); // smax
    double s = (s2 - s1) * sp + s1;

    double t1 = collcor(mchi, z0)*tmin(s, mchi);
    double t2 = tmax();
    double t = -t1 * tp + t1;

    double EstJac = (E2 - E1) * (s2 - s1) * (t2 - t1); // Jacobian to change (Eu, s, t) to unit cube
    // printf("%0.5e\n",tp);

    // double soln = solnFunc(Eu, s, t, uchi, mchi, B, muFe);
    // *HeavisidePhaseSpace(Eu, s, B, s1, s2)
    double int_main = (r * r  * zeta(nE, muFe) * fMB(uchi) * Eu * s * opersdXS(s, t, mchi, oper) * FDfacInitial(Eu, muFe, temp) * FDfacFinal(Eu, s, t, uchi, mchi, B, muFe, temp) / beta_cap(s, mchi) / gamma_cap(s, mchi) / sqrt(B));

    // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\n", Eu, uchi, B, s1);

    return EstJac * int_main;
}

double integrandNoPB_Uappx_T(double r, double uchi, double Ep, double sp, double tp, double muFe, double nE, double B, double mchi, double temp, int oper, double ct)
{

    // No Pauli Blocking Factor

    double E1 = 1.;
    double E2 = EumaxT(muFe, temp, ct);
    double Eu = (E2 - E1) * Ep + E1;

    double s1 = smin(Eu, uchi, B, mchi); // smin
    double s2 = smax(Eu, uchi, B, mchi); // smax
    double s = (s2 - s1) * sp + s1;

    double t1 = tmin(s, mchi);
    double t2 = tmax();
    double t = -t1 * tp + t1;

    double EstJac = (E2 - E1) * (s2 - s1) * (t2 - t1); // Jacobian to change (s,t) to unit cube
    // printf("%0.5e\n",tp);

    double int_main = (r * r * zeta(nE, muFe)  * Eu * s * opersdXS(s, t, mchi, oper) * FDfacInitial( Eu, muFe, temp) / beta_cap(s, mchi) / gamma_cap(s, mchi) / sqrt(B));

    // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\n", Eu, uchi, B, s1);

    return EstJac * int_main;
}

double integrand_Uappx_T(double r, double Ep, double sp, double tp, double muFe, double nE, double B, double mchi, double temp, int oper, double z0, double k0, double ct)
{

    // Withought fMB(u,T)

    double E1 = EuminT(muFe, mchi, temp, k0, ct);
    double E2 = EumaxT(muFe, temp,ct);
    double Eu = (E2 - E1) * Ep + E1;

    double s1 = smin(Eu, 0, B, mchi); // smin
    double s2 = smax(Eu, 0, B, mchi); // smax
    double s = (s2 - s1) * sp + s1;

    double t1 = collcor(mchi, z0)*tmin(s, mchi);
    double t2 = tmax();
    double t = -t1 * tp + t1;

    double EstJac = (E2 - E1) * (s2 - s1) * (t2 - t1); // Jacobian to change (Eu, s, t) to unit cube
    // printf("%0.5e\n",tp);

    // *HeavisidePhaseSpace(Eu, s, B, s1, s2)
    double int_main = (r * r  * zeta(nE, muFe) * Eu * s * opersdXS(s, t, mchi, oper) * FDfacInitial(Eu, muFe, temp) * FDfacFinal(Eu, s, t, 0, mchi, B, muFe, temp) / beta_cap(s, mchi) / gamma_cap(s, mchi) / sqrt(B));

    // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\n", Eu, 0, B, s1);

    return EstJac * int_main;
}

//=============================================================================
// MONTE CAPTURE INTEGRANDS
//=============================================================================

double monteIntegrandNoPB_T(double x[], size_t dim, void *p)
{
    (void)(dim);

    struct crate_params_T *params = (struct crate_params_T *)p;

    double mdm = (params->mchi);
    double T = (params->temp);
    int np = (params->npts);
    int oper = (params->oper);
    double ct = (params->ct);

    double r = x[0]; // r/Rstar
    double u = x[1]; // uchi
    double E = x[2]; // Ep = (Eumax - Eumin)*Eu + Eumin
    double s = x[3]; // sp = (smax - smin)*s + smin
    double t = x[4]; // tp = -tmin*t + tmin

    double muf = muFeinterp(r, np);
    double ne = nEinterp(r, np);
    double B = Binterp(r, np);

    double monteInt = integrandNoPB_T(r, u, E, s, t, muf, ne, B, mdm, T, oper, ct);

    return monteInt;
}

double monteIntegrandFull_T(double x[], size_t dim, void *p)
{
        (void)(dim);

        struct crate_params_T *params = (struct crate_params_T *)p;

        double mdm = (params->mchi);
        double T = (params->temp);
        int np = (params->npts);
        int oper = (params->oper);
        double z0 = (params->z0);
        double k0 = (params->k0);
        double ct = (params->ct);

        double r = x[0]; // r/Rstar
        double u = x[1]; // uchi
        double E = x[2]; // Ep = (Eumax - Eumin)*Eu + Eumin
        double s = x[3]; // sp = (smax - smin)*s + smin
        double t = x[4]; // tp = -tmin*t + tmin

        double muf = muFeinterp(r, np);
        double ne = nEinterp(r, np);
        double B = Binterp(r, np);

        double monteInt = integrand_T(r, u, E, s, t, muf, ne, B, mdm, T, oper, z0, k0, ct);

        return monteInt;
}

double monteIntegrandUappx_T(double x[], size_t dim, void *p)
{
    (void)(dim);

    struct crate_params_T *params = (struct crate_params_T *)p;

    double mdm = (params->mchi);
    double T = (params->temp);
    int np = (params->npts);
    int oper = (params->oper);
    double z0 = (params->z0);
    double k0 = (params->k0);
    double ct = (params->ct);

    double r = x[0]; // r/Rstar
    // double u = x[1]; // uchi
    double E = x[1]; // Ep = (Eumax - Eumin)*Eu + Eumin
    double s = x[2]; // sp = (smax - smin)*s + smin
    double t = x[3]; // tp = -tmin*t + tmin

    double muf = muFeinterp(r, np);
    double ne = nEinterp(r, np);
    double B = Binterp(r, np);

    double monteInt = integrand_Uappx_T(r, E, s, t, muf, ne, B, mdm, T, oper, z0, k0, ct);

    return monteInt;
}

//=====================================================
// INTEGRATORS
//=====================================================

double crateNoPb_T(double mchi, double temp, int oper, int npts, void *cont_vars)
{

    double res, err;

    struct cont_pars *contp  = (struct cont_pars *)cont_vars;

    double z0 = (contp->z0);
    double k0 = (contp->k0);
    double ct = (contp->ct);

    struct crate_params_T params = {mchi, temp, oper, npts, z0, k0, ct};

    double term_max[5] = {1, umax(), 1, 1, 1};
    double term_min[5] = {0, 0, 0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function NoPBfunc = {&monteIntegrandNoPB_T, 5, &params};

    size_t calls = 10000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    {
        gsl_monte_plain_state *s = gsl_monte_plain_alloc(5);
        gsl_monte_plain_integrate(&NoPBfunc, term_min, term_max, 5, calls, rng, s, &res, &err);
        gsl_monte_plain_free(s);

        // display_results("plain", res, err);
    }

    gsl_rng_free(rng);

    return res * Coeff(mchi);
}

double crateFull_T(double mchi, double temp, int oper, int npts, void *cont_vars)
{


    double res, err;

    struct cont_pars *contp  = (struct cont_pars *)cont_vars;

    double z0 = (contp->z0);
    double k0 = (contp->k0);
    double ct = (contp->ct);

    struct crate_params_T params = {mchi, temp, oper, npts, z0, k0, ct};

    double term_max[5] = {1, umax(), 1, 1, 1};
    double term_min[5] = {0, 0, 0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function FullMfunc = {&monteIntegrandFull_T, 5, &params};

    size_t calls = 25000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    if (mchi > 1e1)
    {

        {
            gsl_monte_miser_state *s = gsl_monte_miser_alloc(5);
            gsl_monte_miser_integrate(&FullMfunc, term_min, term_max, 5, calls, rng, s, &res, &err);
            gsl_monte_miser_free(s);

            // display_results("misner", res, err);
        }
    }

    else
    {

        {
            gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(5);

            gsl_monte_vegas_integrate(&FullMfunc, term_min, term_max, 5, calls / 5, rng, s, &res, &err);
            // display_results("vegas warm-up", res, err);
            // printf("converging...\n");

            int MAX_ITERS = 15;
            int iters = 0;

            do
            {
                gsl_monte_vegas_integrate(&FullMfunc, term_min, term_max, 5, calls, rng, s, &res, &err);
                printf("result = % .6e sigma = % .6e "
                        "chisq/dof = %.1e\n",
                        res, err, gsl_monte_vegas_chisq(s));
                printf("at iteration  %d\n", iters);
                iters++;

                if (iters >= MAX_ITERS)
                {
                    printf("Max iterations reached without convergence\n");
                }
            }

            while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5 && iters < MAX_ITERS);

            display_results("vegas final", res, err);

            gsl_monte_vegas_free(s);
        }
    }

    gsl_rng_free(rng);

    return res * Coeff(mchi);
}

double crateUappx_T(double mchi, double temp, int oper, int npts, void *cont_vars)
{
    double res, err;

    struct cont_pars *contp  = (struct cont_pars *)cont_vars;

    double z0 = (contp->z0);
    double k0 = (contp->k0);
    double ct = (contp->ct);

    struct crate_params_T params = {mchi, temp, oper, npts, z0, k0, ct};

    double term_max[4] = {1, 1, 1, 1};
    double term_min[4] = {0, 0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function monteFunc = {&monteIntegrandUappx_T, 4, &params};

    size_t calls = 25000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    if (mchi>1e1){

    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(4);
        gsl_monte_miser_integrate(&monteFunc, term_min, term_max, 4, calls, rng, s, &res, &err);
        gsl_monte_miser_free(s);

        // display_results("misner", res, err);
    }
    }

    else{

    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(4);

        gsl_monte_vegas_integrate(&monteFunc, term_min, term_max, 4, calls/5, rng, s, &res, &err);
        // display_results("vegas warm-up", res, err);
        // printf("converging...\n");

        int MAX_ITERS = 15;
        int iters = 0;

        do
        {
            gsl_monte_vegas_integrate(&monteFunc, term_min, term_max, 4, calls, rng, s, &res, &err);
            printf("result = % .6e sigma = % .6e "
                    "chisq/dof = %.1e\n",
                    res, err, gsl_monte_vegas_chisq(s));
            printf("at iteration  %d\n", iters);
            iters++;

            if (iters >= MAX_ITERS)
            {
                printf("Max iterations reached without convergence\n");
            }
        }

        while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5 && iters < MAX_ITERS);

        // display_results("vegas final", res, err);

        gsl_monte_vegas_free(s);
    }
    }

        gsl_rng_free(rng);

        return res * CoeffUappx(mchi);
    }