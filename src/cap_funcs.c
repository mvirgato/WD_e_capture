#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
// #include <stdint.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>

#include "globals.h"
#include "cap_funcs.h"
#include "cap_interp.h"
// #include "evap_funcs.h"

//=============================================================================
// MISC. FUNCTIONS
//=============================================================================

void display_results (char *title, double result, double error){
    printf ("%s ==================\n", title);
    printf ("result = % 0.5e\n", result);
    printf ("sigma  = % 0.5e\n\n", error);
}

//=============================================================================
// SYSTEM SELECT
//=============================================================================
void system_setup(double mT, char *WDname, char *eos_type, char *elem, char *vstar, int *np)
{
    // mT sets targer mass, enter EoS file to read in

    char *EOSname;

    if (eos_type == "r"){
        EOSname = malloc(strlen("./EoSs/FMT_") + strlen(WDname) + strlen("e6m_" )+ strlen(elem) + strlen(".dat") + 1);
        sprintf(EOSname, "./EoSs/FMT_%se6m_%s.dat", WDname, elem);
    }
    else if (eos_type == "M"){
        EOSname = malloc(strlen("./EoSs/FMT_") + strlen(WDname) + strlen("_Msun_")+ strlen(elem) + strlen(".dat") + 1);
        sprintf(EOSname, "./EoSs/FMT_%s_Msun_%s.dat", WDname, elem);
    }

    char *mstarname = malloc(strlen("./mstar_profs/mstar_e_WD_") + strlen(vstar) + strlen(".dat") + 1);
    sprintf(mstarname, "./mstar_profs/mstar_e_WD_%s.dat", vstar);

    read_profs(EOSname, np);
    read_mstar(mstarname); // need a way to change for different mass stars
    mt = mT;
    free(EOSname);
    free(mstarname);
}

//=====================================================
// CAPTURE FUNCTIONS
//=====================================================

double fMB(double uchi){
    // Maxwell Boltzmann distribution  x uchi without constants

    return (exp(-3. * (uchi - vs * kmTOm / cspeed) * (uchi - vs * kmTOm / cspeed) / (2. * vd * vd * kmTOm * kmTOm / cspeed / cspeed)) - exp(-3. * (uchi + vs * kmTOm / cspeed) * (uchi + vs * kmTOm / cspeed) / (2. * vd * vd * kmTOm * kmTOm / cspeed / cspeed)));
}

double beta_cap(double s, double mchi){
    // beta in cross sections
    return (s - mt*mt - mchi*mchi);
}

double gamma_cap(double s, double mchi){
    // gamma in cross sections
    double b = beta_cap(s, mchi);

    return sqrt(b*b - 4*mchi*mchi*mt*mt);
}

double solnFunc(double Eu, double s, double t, double uchi, double mchi, double B, double muFe){

    double E1 = mchi*(1 + uchi*uchi/2.)/sqrt(B);
    double E2 = mt*Eu;

    double mchi2 = mchi*mchi;
    double mchi4 = mchi2*mchi2;
    double mchi6 = mchi2*mchi2*mchi2;

    double mt2 = mt*mt;
    double mt4 = mt2*mt2;
    double mt6 = mt2*mt2*mt2;

    double E12 = E1*E1;
    double E22 = E2*E2;

    double soln;

    if ( s > mt2 + mchi2 + 2.*E1*E2){
        soln = (-(E2 * mchi4 * t) + E2 * mt4 * t - 2 * E2 * mt2 * s * t + E2 * s*s * t + 2 * E12 * E2 * (mt4 + (mchi2 - s) * (mchi2 - s - t) + mt2 * (-2 * mchi2 - 2 * s + t)) + sqrt(pow(2 * E1 * E2 + mchi2 + mt2 - s, 2) * (4 * E22 * mchi2 + mchi4 + 4 * E12 * mt2 + mt4 + 4 * E1 * E2 * (mchi2 + mt2 - s) -2 * mchi2 * s + s*s - 2 * mt2 * (mchi2 + s)) *t * (mchi4 + mt4 - 2 * mchi2 * s - 2 * mt2 * (mchi2 + s) + s * (s + t))) + E1 * (pow(mchi, 6) + pow(mt, 6) + mt4 * (-mchi2 - 3 * s + t) - mchi4 * (3 * s + t) - mt2 * (mchi4 + 2 * mchi2 * s - 3 * s*s - 2 * E22 * t) - s * (s*s + 2 * E22 * t + s * t) + mchi2 * (3 * s*s - 2 * E22 * t + 2 * s * t))) /((2 * E1 * E2 + mchi2 + mt2 - s) * (mt4 + pow(mchi2 - s, 2) - 2 * mt2 * (mchi2 + s)));
    }
    else{
        soln = (-(E2 * mchi4 * t) + E2 * mt4 * t - 2 * E2 * mt2 * s * t + E2 * s*s * t + 2 * E12 * E2 * (mt4 + (mchi2 - s) * (mchi2 - s - t) + mt2 * (-2 * mchi2 - 2 * s + t)) - sqrt(pow(2 * E1 * E2 + mchi2 + mt2 - s, 2) * (4 * E22 * mchi2 + mchi4 + 4 * E12 * mt2 + mt4 + 4 * E1 * E2 * (mchi2 + mt2 - s) - 2 * mchi2 * s + s*s - 2 * mt2 * (mchi2 + s)) * t * (mchi4 + mt4 - 2 * mchi2 * s - 2 * mt2 * (mchi2 + s) + s * (s + t))) + E1 * (pow(mchi, 6) + pow(mt, 6) + mt4 * (-mchi2 - 3 * s + t) - mchi4 * (3 * s + t) - mt2 * (mchi4 + 2 * mchi2 * s - 3 * s*s - 2 * E22 * t) - s * (s*s + 2 * E22 * t + s * t) + mchi2 * (3 * s*s - 2 * E22 * t + 2 * s * t))) / ((2 * E1 * E2 + mchi2 + mt2 - s) * (mt4 + pow(mchi2 - s, 2) - 2 * mt2 * (mchi2 + s)));
    }

    return soln;
}

double HeavisidePB(double Eu, double s, double t, double uchi, double mchi, double B, double muFe, double soln){

    double E1 = mchi*(1 + uchi*uchi/2.)/sqrt(B);
    double E2 = mt*Eu;

    double mchi2 = mchi*mchi;
    double mchi4 = mchi2*mchi2;
    double mchi6 = mchi2*mchi2*mchi2;

    double mt2 = mt*mt;
    double mt4 = mt2*mt2;
    double mt6 = mt2*mt2*mt2;

    double E12 = E1*E1;
    double E22 = E2*E2;

    double EF = E2 + E1 - soln - mt - muFe;

    if (isnan(EF) == 0){
        if (EF>0){
            return 1.;
        }
        else{
            return 0.;
        }
    }
    else{
        return 0.;
    }
}

double HeavisideUchi(double Eu, double s, double t, double uchi, double mchi, double B, double muFe, double soln){

    double E1 = mchi*(1 + uchi*uchi/2.)/sqrt(B);
    double E2 = mt*Eu;

    double mchi2 = mchi*mchi;
    double mchi4 = mchi2*mchi2;
    double mchi6 = mchi2*mchi2*mchi2;

    double mt2 = mt*mt;
    double mt4 = mt2*mt2;
    double mt6 = mt2*mt2*mt2;

    double E12 = E1*E1;
    double E22 = E2*E2;

    double finite_q0 = E1 - soln - 0.5*mchi*uchi*uchi/sqrt(B);

    if (isnan(finite_q0) == 0){
        if (finite_q0>0){
            return 1.;
        }
        else{
            return 0.;
        }
    }
    else{
        return 0.;
    }
}

double HeavisidePhaseSpace(double Eu, double s, double B, double _smin, double _smax){
    double cos = (-2.*s + _smax + _smin)/( _smax - _smin);

    double condition = Eu * sqrt( (1 - B)/(Eu*Eu - 1) ) - cos;

    if (condition > 0){
        return 1.;
    }
    else {
        return 0.;
    }
}



double nfree(double muFe){
    // free number density
    double num = (muFe * (2. * mt + muFe)) * sqrt(muFe * (2. * mt + muFe));
    double denom = 3. * pi * pi ;

    return num/denom;
}

double zeta(double nE, double muFe){
    // correction to free number density

    double nf0 = nfree(muFe0);
    double zeta0 =  (n0* 3. * pi * pi * inveVTOm * mTOpm * inveVTOm * mTOpm * inveVTOm * mTOpm) / (muFe0 * (2. * mt + muFe0) * sqrt(muFe0 * (2. * mt + muFe0)));
    return (nE*zeta0*nf0)/(n0*nfree(muFe));
}

double Coeff(double mchi){
    // numerical factors

    return (2. * rhoDM * sqrt(3. / 2. / pi) * (RS * RS * RS * mTOinveV * mTOinveV * mTOinveV) * mt * mt) / ((vs * vd * kmTOm * kmTOm / cspeed / cspeed) * (cmTOm * cmTOm * cmTOm * mTOinveV * mTOinveV * mTOinveV) * pi * mchi * mchi * inveVTOs);
}


//=====================================================
// INTEGRATION LIMITS
//=====================================================

double collcor(double mchi, double z0){
    return 1. - fmax(1. - z0 * mchi / mt, 0);
}

double tmin(double s, double mchi){
    return -(gamma_cap(s, mchi) * gamma_cap(s, mchi))/s;
}

double tmax(){
    return 0.;
}

double smin(double Eu, double uchi, double B, double mchi){
    return (mt * mt + mchi * mchi + (2. * Eu * mt * mchi * (1. + (uchi * uchi) / 2.)) / sqrt(B) - 2. * sqrt(-1. + (1. + uchi * uchi / 2) * (1. + uchi * uchi / 2) / B)*sqrt(-1. + Eu*Eu)*mt*mchi);
}

double smax(double Eu, double uchi, double B, double mchi){
    return (mt * mt + mchi * mchi + (2. * Eu * mt * mchi * (1. + (uchi * uchi) / 2.)) / sqrt(B) + 2. * sqrt(-1. + (1. + uchi * uchi / 2) * (1. + uchi * uchi / 2) / B) * sqrt(-1. + Eu*Eu) * mt * mchi);
}

double Eumin(double muFe, double mchi, double k0){
    return fmax(1., 1+(muFe/mt)*(1 - k0*mchi/mt) );
}

double Eumax(double muFe){
    return (1. + muFe/mt);
}

double umax(){
    return 10.*vd*kmTOm/cspeed;
}

double rmin(double mchi){
    double lm = log10(mchi);
    return fmax(0., 0.9 - (pow(10, (-1.35*lm)) / 1e9 + pow(10,(1.38*lm)) / pow(10,2.7)));
    // return fmax(0., 0.8 - (pow(10,(1.38*lm)) / pow(10,1)));
    // return 0.;

}

double rmax(){
    return 1.;
}

//=====================================================
// CAPTURE INTEGRANDS
//=====================================================

double integrandNoPB(double r, double uchi, double Ep, double sp, double tp, double muFe, double nE, double B, double mchi, int oper){

    // No Pauli Blocking Factor

    double E1 = 1.;
    double E2 = Eumax(muFe);
    double Eu = (E2 - E1)*Ep + E1;

    double s1 = smin(Eu, uchi, B, mchi); // smin
    double s2 = smax(Eu, uchi, B, mchi); // smax
    double s  = (s2- s1)*sp + s1;

    double t1 = tmin(s, mchi);
    double t2 = tmax();
    double t  = -t1*tp + t1; 

    double EstJac = (E2 - E1) * (s2 - s1) * (t2 - t1); // Jacobian to change (s,t) to unit cube
    // printf("%0.5e\n",tp);

    double soln = solnFunc(Eu, s, t, uchi, mchi, B, muFe);
  
    double int_main = (r * r * zeta(nE, muFe) * fMB(uchi) * Eu * s * opersdXS(s, t, mchi, oper) * HeavisideUchi(Eu, s, t, uchi, mchi, B, muFe, soln)/ beta_cap(s, mchi) / gamma_cap(s, mchi) / sqrt(B));

    // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\n", Eu, uchi, B, s1);

    return EstJac*int_main;
}

double integrandFull(double r, double uchi, double Ep, double sp, double tp, double muFe, double nE, double B, double mchi, int oper, double z0, double k0)
{
    // full integrand

    double E1 = Eumin(muFe, mchi, k0);
    double E2 = Eumax(muFe);
    double Eu = (E2 - E1) * Ep + E1;

    double s1 = smin(Eu, uchi, B, mchi); // smin
    double s2 = smax(Eu, uchi, B, mchi); // smax
    double s = (s2 - s1) * sp + s1;

    double t1 = collcor(mchi, z0)*tmin(s, mchi);
    double t2 = tmax();
    double t = -t1 * tp + t1;

    double EstJac = (E2 - E1) * (s2 - s1) * (t2 - t1); // Jacobian to change (s,t) to unit cube
    // printf("%0.5e\n",tp);

    double soln = solnFunc(Eu, s, t, uchi, mchi, B, muFe);
    //  * HeavisidePhaseSpace(Eu, s, B, s1, s2)
    double int_main = (r * r * zeta(nE, muFe) * fMB(uchi) * Eu * s * opersdXS(s, t, mchi, oper)*HeavisidePB(Eu, s, t, uchi, mchi, B, muFe, soln) * HeavisideUchi(Eu, s, t, uchi, mchi, B, muFe, soln) / beta_cap(s, mchi) / gamma_cap(s, mchi) / sqrt(B));

    // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\n", Eu, uchi, B, s1);

    // printf("%f\n", r);
    return EstJac * int_main;
}






//=============================================================================
// MONTE CAPTURE INTEGRANTS
//=============================================================================


double monteIntegrandNoPB(double x[], size_t dim, void * p){
    (void)(dim);

    struct crate_params *params = (struct crate_params *)p;

    double mdm   = (params->mchi);
    int    np    = (params->npts);
    int    oper  = (params->oper);

    double r = x[0]; // r/Rstar
    double u = x[1]; // uchi
    double E = x[2]; // Ep = (Eumax - Eumin)*Eu + Eumin
    double s = x[3]; // sp = (smax - smin)*s + smin
    double t = x[4]; // tp = -tmin*t + tmin

    double muf = muFeinterp(r, np);
    double ne  = nEinterp(r, np);
    double B   = Binterp(r, np);
    
    double monteInt = integrandNoPB(r, u, E, s, t, muf, ne, B, mdm, oper);

    return monteInt;

}

double monteIntegrandFull(double x[], size_t dim, void *p)
{
    (void)(dim);

    struct crate_params *params = (struct crate_params *)p;

    double mdm = (params->mchi);
    int np = (params->npts);
    int oper = (params->oper);
    double z0 = (params->z0);
    double k0 = (params->k0);

    double r = x[0]; // r/Rstar
    double u = x[1]; // uchi
    double E = x[2]; // Ep = (Eumax - Eumin)*Eu + Eumin
    double s = x[3]; // sp = (smax - smin)*s + smin
    double t = x[4]; // tp = -tmin*t + tmin
    

    double muf = muFeinterp(r, np);
    double ne = nEinterp(r, np);
    double B = Binterp(r, np);

    double monteInt = integrandFull(r, u, E, s, t, muf, ne, B, mdm, oper, z0, k0);

    return monteInt;
}



//=====================================================
// INTEGRATORS
//=====================================================

double crateNoPb(double mchi, int oper, int npts){

    double res, err;

    struct crate_params params = {mchi, oper, npts};

    double term_max[5] = {1, umax(), 1, 1, 1};
    double term_min[5] = {0, 0, 0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function NoPBfunc = {&monteIntegrandNoPB, 5, &params};

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

    return 0.5*res*Coeff(mchi);
}


double crateFull(double mchi, int oper, int npts, void *cont_vars)
{


    double res, err;

    struct cont_pars *contp  = (struct cont_pars *)cont_vars;

    double z0 = (contp->z0);
    double k0 = (contp->k0);

    struct crate_params params = {mchi, oper, npts, z0, k0};

    double term_max[5] = {1, umax(), 1, 1, 1};
    double term_min[5] = {0, 0, 0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function FullMfunc = {&monteIntegrandFull, 5, &params};

    size_t calls = 20000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    // {
    //     gsl_monte_plain_state *s = gsl_monte_plain_alloc(5);
    //     gsl_monte_plain_integrate(&FullMfunc, term_min, term_max, 5, calls, rng, s, &res, &err);
    //     gsl_monte_plain_free(s);

    //     // display_results("plain", res, err);
    // }

    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(5);
        gsl_monte_miser_integrate(&FullMfunc, term_min, term_max, 5, calls, rng, s, &res, &err);
        gsl_monte_miser_free(s);

        // display_results("misner", res, err);
    }

    // {

    //     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(5);

    //     gsl_monte_vegas_integrate(&FullMfunc, term_min, term_max, 5, calls/5, rng, s, &res, &err);
    //     // display_results("vegas warm-up", res, err);
    //     // printf("converging...\n");

        // int MAX_ITERS = 15;
        // int iters = 0;
        

    //     do
    //     {
    //         gsl_monte_vegas_integrate(&FullMfunc, term_min, term_max, 5, calls, rng, s, &res, &err);
    //         printf("result = % .6e sigma = % .6e "
    //                "chisq/dof = %.1e\n",
    //                res, err, gsl_monte_vegas_chisq(s));
    //         printf("at iteration  %d\n", iters);
    //         iters++; 

    //         if (iters >= MAX_ITERS)
    //         {
    //             printf("Max iterations reached without convergence\n");
    //         }
    //     }

    //     while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5 && iters < MAX_ITERS);

    //     display_results("vegas final", res, err);

    //     gsl_monte_vegas_free(s);
    // }

    gsl_rng_free(rng);

    return res * Coeff(mchi);
}



//=====================================================
// CROSS SECTIONS
//=====================================================

double dXS_const(double s, double t, double mchi){
    // constant differential cross section

    (void)(s); (void)(t); (void)(mchi);

    return 1e-39 * mTOinveV*mTOinveV*cmTOm*cmTOm/2.;
}

double opersdXS(double s, double t, double mchi, int num){
    // enter num = 11 for constant cross section
    
    double prefac = beta_cap(s, mchi) / 16. / pi / (2. * s * beta_cap(s, mchi) -  gamma_cap(s, mchi) * gamma_cap(s, mchi));

    double mu  = mchi/mt;
    double mu2 = mu*mu;
    double mu4 = mu*mu*mu*mu;
    
    double mchi2  = mchi*mchi;
    double mchi4  = mchi*mchi*mchi*mchi;

    double mt2 = mt*mt;
    double mt4 = mt * mt * mt * mt;

    double M2;

    if (num == 1){
        M2 = (-4 * mchi2 + t) * (-4 * mt2 + t);
        // printf("in num %d\n", num);
    }
    else if (num == 2){
        M2 = t * (-4 * mt2 + t);
        // printf("in num %d\n", num);
    }
    else if (num == 3){
        M2 = t * (-4 * mchi2 + t);
        // printf("in num %d\n", num);
    }
    else if (num == 4){
        M2 =  t * t;
        // printf("in num %d\n", num);
    }
    else if (num == 5){
        M2 = 2 * (2 * mchi4 + 2 * mt4 + 4 * mchi2 * (mt2 - s) - 4 * mt2 * s + 2 * s*s + 2 * s * t + t*t);
    }
    else if (num == 6){
        M2 = 2 * (2 * mchi4 + 2 * mt4 - 4 * mt2 * s + 2 * s*s + 2 * s * t + t*t - 4 * mchi2 * (mt2 + s + t));
    }
    else if (num == 7){
        M2 = 2 * (2 * mchi4 + 2 * mt4 + 2 * s*s - 4 * mchi2 * (mt2 + s) + 2 * s * t + t*t - 4 * mt2 * (s + t));
    }
    else if (num == 8){
        M2 = 2 * (2 * mchi4 + 2 * mt4 + 2 * s*s + 4 * mchi2 * (5 * mt2 - s - t) + 2 * s * t + t*t - 4 * mt2 * (s + t));
    }
    else if (num == 9){
        M2 = 8 * (4 * mchi4 + 4 * mt4 + 2 * mchi2 * (8 * mt2 - 4 * s - t) + pow(2 * s + t, 2) - 2 * mt2 * (4 * s + t));
    }
    else if (num == 10){
        M2 = 8 * (4 * mchi4 + 4 * mt4 + pow(2 * s + t, 2) - 2 * mt2 * (4 * s + t) - 2 * mchi2 * (4 * mt2 + 4 * s + t));
    }
    else if (num == 11){
        M2 = 1e-39 * mTOinveV*mTOinveV*cmTOm*cmTOm/2; // const. XS
        // M2 = t*t;
    }

    // double M2l[11] = {

    //     (4. * mchi * mchi - t) * (4. * mchi * mchi - mu * mu * t) / mu / mu,
    //     t * (mu * mu * t - 4. * mchi * mchi) / mu / mu,
    //     t * (t - 4. * mchi * mchi),
    //     t * t,
    //     2. * (2. * (mu * mu + 1.) * (mu * mu + 1.) * mchi4 - 4. * (mu2 + 1.) * mu2 * s * mchi2 + mu4 * (2. * s * s + 2. * s * t + t * t)) / mu4,
    //     2. * (2. * (mu2 - 1) * (mu2 - 1.) * mchi4 - 4. * mu2 * mchi2 * (mu2 * s + s + mu2 * t) + mu4 * (2. * s * s + 2. * s * t + t * t)) / mu4,
    //     2. * (2. * (mu2 - 1) * (mu2 - 1.) * mchi4 - 4. * mu2 * mchi2 * (mu2 * s + s + t) + mu4 * (2. * s * s + 2. * s * t + t * t)) / mu4,
    //     2. * (2. * (mu4 + 10. * mu2 + 1.) * mchi4 - 4.*(mu2 + 1.) * mu2 * mchi2 * (s + t) + mu4 * (2. * s * s + 2. * s * t + t * t)) / mu4,
    //     8. * (4. * (mu4 + 4. * mu2 + 1.) * mchi4 - 2. * (mu2 + 1.) * mu2 * mchi2 * (4. * s + t) + mu4 * (2. * s + t) * (2. * s + t)) / mu4,
    //     8. * (4. * (mu2 - 1.) * (mu2 - 1.) * mchi4 - 2. * (mu2 + 1.) * mu2 * mchi2 * (4. * s + t) + mu4 * (2. * s + t) * (2. * s + t)) / mu4,
    //     1e-39 * mTOinveV*mTOinveV*cmTOm*cmTOm/2. // const. XS
    // };

    double Ci[11] = {2 * mt * mt / vev / vev, 2 * mt * mt / vev / vev, 2 * mt * mt / vev / vev, 2 * mt * mt / vev / vev, 1, 1, 1, 1, 1, 1, 1}; // Operator couplings

    if (num<=10 && num>=1)
    {
        return prefac * Ci[num-1] * M2;
    }
    else if(num == 11){
        return prefac * M2;
    }
    else
    {
        printf("Enter number from 1-10\nfor corresponding operator\nor 11 for const. cross section");
        exit(0);
        return 0;
    }

    
}
