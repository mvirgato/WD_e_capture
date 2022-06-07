#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_erf.h>

#include "globals.h"
#include "cap_funcs.h"
#include "cap_interp.h"
#include "evap_funcs.h"

//=============================================================================
// MISC FUNCTIONS
//=============================================================================

//=============================================================================
// EVAPORATION FUNCTIONS
//=============================================================================

double radChiSqr(double T, double mchi){
    // r^2 = 3 T / 2 G pi rho_c n_chi
    //T in K, mchi in MeV
    // output in m^3

    return 3. * T * kBMeV*cspeed*cspeed/ 2./pi/rhoc/mchi/GNewt;
}

double niso(double rad, double T, double mchi){
    // normalisation is in the prefactor definition

    double rchi2 = radChiSqr(T, mchi);
    double norm = pi * rchi2 * (-2. * exp(-RS * RS / rchi2) + sqrt(pi) * sqrt(rchi2) * erf(RS / sqrt(rchi2)));

    return  exp(-rad*rad*RS*RS/rchi2)/norm;
}

double fMBevap(double vel, double vesc, double T, double mchi){
    // T in K

    double velchi2 = kBMeV*T*2./mchi;

    double velMB = exp(-vel*vel/velchi2);

    double norm = pi * velchi2 * cspeed * cspeed * cspeed * (sqrt(pi) * sqrt(velchi2) * erf(vesc / sqrt(velchi2)) - 2. * vesc * exp(-vesc * vesc / velchi2) );

    return velMB;
}


double heavisideprod(double w, double v, double t, double s){

    if ((t + s - v)>0 && (w - fabs(t-s))>0){
       return 1.; 
    }
    else{
        return 0.;
    }
}

double prefactorEvap(double mchi, double T){
    double mup = (1. + mchi/mt)/2.;
    double mup4 = mup*mup*mup*mup;
    double k2 = mt/2./T;
    double k3 = k2*sqrt(k2);

    return (512. * pi * sqrt(pi) * mup4 * k3 * RS * RS * RS * cspeed * cspeed) / (pmTOm * pmTOm * pmTOm );
}

//=============================================================================
// INTEGRANDS
//=============================================================================

double evapIntegrand(double r, double wp, double vp, double tp, double sp, double mchi, double vesc, double T, double ne) {
     
    double mu = mchi/mt;
    double mup = (1. + mu)/2.;
    double k2 = mt/2./T;

    double w = vesc*wp;
    double v = vesc + (1. - vp) / vp;
    double t = (1. - tp) / tp;
    double s = (1. - sp) / sp;

    double Jac = vesc/vp/vp/tp/tp/sp/sp;

    // double v = ((1. + 0.3)*vesc - vesc)*vp - vesc;
    // double x = ((1. + 0.2)*v - v)*sp - v;
    // double y = (2.*w)*tp - w;
    // double s = (x + y)/2.;
    // double t = (x - y)/2.;

    // double Jac = ((1 + 0.3)*vesc - vesc)*vesc*((1. + 0.2)*v - v)*(2.*w)/2.;

    double velec2 = 2. * mu * mup * t * t + 2. * mup * s * s - mu * w * w;


    double intmain = r * r * niso(r, T, mchi) * fMBevap(w, vesc, T, mchi) * ne * v * t * w * exp(-k2 * velec2) * heavisideprod(w, v, t, s);

    if (intmain != 0){
        printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n", s, t, w, niso(r, T, mchi), intmain);
    }
    
    return Jac*intmain;

}

double reg1evapIntegrand(double r, double wp, double vp, double tp, double sp, double mchi, double vesc, double T, double ne)
{

    double mu = mchi / mt;
    double mup = (1. + mu) / 2.;
    double k2 = mt / 2. / T;

    double w = vesc * wp;
    double v = vesc + (1. - vp) / vp;

    double t1 = (v - w)/2.;
    double t2 = (v + w)/2.;
    double t = (t2 - t1)*tp + t1;

    double s1 = v - t;
    double s2 = w + t;
    double s = (s2 - s1)*sp + s1;

    double velec2 = 2. * mu * mup * t * t + 2. * mup * s * s - mu * w * w;

    double Jac = vesc * ( t2 - t1) * (s2 - s1)/vp/vp;


    double intmain = r * r * niso(r, T, mchi) * fMBevap(w, vesc, T, mchi) * ne * v * t * exp(-k2 * velec2) / w;

    // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n", r, w, v, t, s, ne);
    return Jac * intmain;
}

//=============================================================================
// MONTE INTEGRANTS
//=============================================================================

double monteEvapIntegrand(double x[], size_t dim, void * p){
    (void)(dim);

    struct evap_params *params = (struct evap_params *)p;

    double mdm  = (params->mchi);
    double T    = (params->T);
    int    np   = (params->npts);

    double rp = x[0];
    double wp = x[1];
    double vp = x[2];
    double tp = x[3];
    double sp = x[4];

    double ne   = nEinterp(rp, np);
    double vesc = sqrt(1. - Binterp(rp, np));

    double monteInt = evapIntegrand(rp, wp, vp, tp, sp, mdm, vesc, T, ne);

    // if (isnan(monteInt) == 1 ){
    //     printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n", rp, wp, vp, tp, sp, monteInt);
    // }
    return monteInt;

}

double monteReg1Integrand(double x[], size_t dim, void *p)
{
    (void)(dim);

    struct evap_params *params = (struct evap_params *)p;

    double mdm = params->mchi;
    double T = params->T;
    int np = params->npts;

    double rp = x[0];
    double wp = x[1];
    double vp = x[2];
    double tp = x[3];
    double sp = x[4];

    double ne = nEinterp(rp, np);
    double vesc = sqrt(1. - Binterp(rp, np));

    double monteInt = reg1evapIntegrand(rp, wp, vp, tp, sp, mdm, vesc, T, ne);

    // if (monteInt == 0 ){
    //     printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n", rp, wp, vp, tp, sp, monteInt);
    // }
    return monteInt;
}

//=====================================================
// INTEGRATORS
//=====================================================

double evaprate(double mchi, double T, int npts){

    double res, err;

    struct evap_params params = { mchi, T, npts };

    double term_max[5] = {1,1,1,1,1};
    double term_min[5] = {0,0,0,0,0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function evapint = {&monteEvapIntegrand, 5, &params};

    size_t calls = 20000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    // {
    //     gsl_monte_plain_state *s = gsl_monte_plain_alloc(5);
    //     gsl_monte_plain_integrate(&evapint, term_min, term_max, 5, calls, rng, s, &res, &err);
    //     gsl_monte_plain_free(s);

    //     display_results("plain", res, err);
    // }

    // {
    //     gsl_monte_miser_state *s = gsl_monte_miser_alloc(5);
    //     gsl_monte_miser_integrate(&evapint, term_min, term_max, 5, calls, rng, s, &res, &err);
    //     gsl_monte_miser_free(s);

    //     display_results("misner", res, err);
    // }

    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(5);

        gsl_monte_vegas_integrate(&evapint, term_min, term_max, 5, calls/10, rng, s, &res, &err);
        display_results("vegas warm-up", res, err);
        printf("converging...\n\n");

        do
        {
            gsl_monte_vegas_integrate(&evapint, term_min, term_max, 5, calls, rng, s, &res, &err);
            printf("result = % .6e sigma = % .6e "
                   "chisq/dof = %.1e\n",
                   res, err, gsl_monte_vegas_chisq(s));
        }

        while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);

        display_results("vegas final", res, err);

        gsl_monte_vegas_free(s);
    }


    gsl_rng_free(rng);

    return res*prefactorEvap(mchi, T)*dXS_const(1,1,mchi);
}

double reg1evap(double mchi, double T, int npts)
{

    double res, err;

    struct evap_params params = {mchi, T, npts};

    double term_max[5] = {1, 1, 1, 1, 1};
    double term_min[5] = {0, 0, 0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function reg1int = {&monteReg1Integrand, 5, &params};

    size_t calls = 10000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    {
        gsl_monte_plain_state *s = gsl_monte_plain_alloc(5);
        gsl_monte_plain_integrate(&reg1int, term_min, term_max, 5, calls, rng, s, &res, &err);
        gsl_monte_plain_free(s);

        display_results("plain", res, err);
    }

    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(5);
        gsl_monte_miser_integrate(&reg1int, term_min, term_max, 5, calls, rng, s, &res, &err);
        gsl_monte_miser_free(s);

        display_results("misner", res, err);
    }

    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(5);

        gsl_monte_vegas_integrate(&reg1int, term_min, term_max, 5, calls/10, rng, s, &res, &err);
        display_results("vegas warm-up", res, err);
        printf("converging...\n\n");

        do
        {
            gsl_monte_vegas_integrate(&reg1int, term_min, term_max, 5, calls, rng, s, &res, &err);
            printf("result = % .6e sigma = % .6e "
                   "chisq/dof = %.1e\n",
                   res, err, gsl_monte_vegas_chisq(s));
        }

        while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);

        display_results("vegas final", res, err);

        gsl_monte_vegas_free(s);
    }

    gsl_rng_free(rng);

    return res * prefactorEvap(mchi, T);
}
