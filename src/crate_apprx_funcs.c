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
#include "multiscatter.h"

//=====================================================
// GENERAL FUNCTIONS
//=====================================================

double CoeffUappx(double mchi)
{
    return (2. * rhoDM * erf(sqrt(3. / 2.) * vs / vd) * mt * mt * (RS * RS * RS * mTOinveV * mTOinveV * mTOinveV)) / (pi * (vs * kmTOm / cspeed) * mchi * mchi * cmTOm * cmTOm * cmTOm * mTOinveV * mTOinveV * mTOinveV * inveVTOs);
}

double coeff_high(double mchi)
{
    return pi*(rhoDM/cmTOm/cmTOm/cmTOm/mTOinveV/mTOinveV/mTOinveV)*erf(sqrt(3./2.)*vs/vd)*(RS*RS*RS)/(2.*(vs*kmTOm / cspeed)*mchi*pmTOm*pmTOm*pmTOm*inveVTOs);
}

// double coeff_high(double mchi)
// {
//     return pi*(rhoDM/cmTOm/cmTOm/cmTOm/mTOinveV/mTOinveV/mTOinveV)*erf(sqrt(3./2.)*vs/vd)*(RS*RS*RS)/(2.*(vs*kmTOm / cspeed)*mchi*mchi*pmTOm*pmTOm*pmTOm*inveVTOs);
// }

//=====================================================
// CROSS SECTIONS FOR HIGH MASS
//=====================================================

// double int_rate_high(double mchi, double r, double ne, double B, double muf, int oper){
double int_rate_high(double mchi, double B, double muf, int oper){

    double mu = mchi / mt;
    double mu2 = mu * mu;
    double mu4 = mu * mu * mu * mu;

    double mchi2 = mchi * mchi;
    double mchi4 = mchi * mchi * mchi * mchi;

    double mt2 = mt * mt;
    double mt4 = mt * mt * mt * mt;

    double muf2 = muf*muf;

    double s0 = mchi2 + 2.*(mt+muf)*mchi/sqrt(B);
    double s02 = s0*s0;
    // double s02 = mchi4 + 4.*(mt+muf)*mchi2*mchi/sqrt(B);

    double beta = s0 - (mt2 + mchi2);
    // double beta = 2.*(mt + muf)*mchi/sqrt(B);

    double gamma2 = beta * beta - 4. * mt2 * mchi2;
    // double gamma2 = -4.*mt2*mchi2 + (4.*mt2*mchi2)/B + (8.*mt*muf*mchi2)/B + (4.*muf2*mchi2)/B;


    double g_list[3];

    if (oper==1){
        g_list[0] = 16.*mchi2*mt2;
        // g_list[1] = (-4.*mchi2 - 4.*mt2);
        g_list[1] = (-4.*mchi2);
        g_list[2] = 1.;
    }
    else if (oper == 2){
        g_list[0] = 0.;
        g_list[1] = -4.*mt2;
        g_list[2] = 1.;
    }
    else if (oper == 3){
        g_list[0] = 0.;
        g_list[1] = -4.*mchi2;
        g_list[2] = 1.;
    }
    else if (oper == 4){
        g_list[0] = 0.;
        g_list[1] = 0.;
        g_list[2] = 1.;
    }
    else if (oper == 5){
        // g_list[0] = 2.*(2.*mchi4 + 2.*mt4 + 4.*mchi2*(mt2 - s0) - 4*mt2*s0 + 2.*s02);
        // g_list[1] = 4.*s0;
        // g_list[2] = 2.;

        g_list[0] = (16. * mchi2 * mt2) / B + (32. * mchi2 * mt * muf) / B + (16. * mchi2 * muf*muf) / B;
        g_list[1] = 4. * mchi2 + (8. * mchi * mt) / sqrt(B) + (8. * mchi * muf) / sqrt(B);
        g_list[2] = 2.;
    }
    else if (oper == 6){
        // g_list[0] = 2.*(2.*mchi4 - 4.*mchi2*mt2 + 2.*mt4 - 4.*mchi2*s0 - 4.*mt2*s0 + 2.*s02);
        // g_list[1] = 2.*(-4.*mchi2 + 2.*s0);
        // g_list[2] = 2.;

        g_list[0] = (-16.*(-1. + B)*mchi2*mt2)/B + (32.*mchi2*mt*muf)/B + (16.*mchi2*muf2)/B;
        g_list[1] =-4.*mchi2 + (8.*mchi*mt)/sqrt(B) + (8.*mchi*muf)/sqrt(B);
        g_list[2] = 2.;
    }
    else if (oper == 7){
        // g_list[0] = 2.*(2.*mchi4 + 2.*mt4 - 4.*mt2*s0 + 2.*s02 - 4.*mchi2*(mt2 + s0)); 
        // g_list[1] = 2.*(-4.*mt2 + 2.*s0); 
        // g_list[2] = 2.;

        g_list[0] = (-16.*(-1. + B)*mchi2*mt2)/B + (32.*mchi2*mt*muf)/B + (16.*mchi2*muf2)/B;
        g_list[1] = 4.*mchi2 + (8.*mchi*mt)/sqrt(B) + (8.*mchi*muf)/sqrt(B);
        g_list[2] = 2.;
    }
    else if (oper == 8){
        // g_list[0] = 4.*mchi4 + 40.*mchi2*mt2 + 4.*mt4 - 8.*mchi2*s0 - 8.*mt2*s0 + 4.*s02; 
        // g_list[1] = -8.*mchi2 - 8.*mt2 + 4.*s0; 
        // g_list[2] = 2.;

        g_list[0] = (16.*(1. + 2.*B)*mchi2*mt2)/B + (32.*mchi2*mt*muf)/B + (16.*mchi2*muf2)/B;
        g_list[1] = -4.*mchi2 + (8.*mchi*mt)/sqrt(B) + (8.*mchi*muf)/sqrt(B);
        g_list[2] = 2.;
    }
    else if (oper == 9){
        g_list[0] = 8*(4*mchi4 + 16*mchi2*mt2 + 4*mt4 - 8*mchi2*s0 - 8*mt2*s0 + 4*s02); 
        g_list[1] = 8*(-2*mchi2 - 2*mt2 + 4*s0); 
        g_list[2] = 8.;

        // g_list[0] = (16.*(1. + 2.*B)*mchi2*mt2)/B + (32.*mchi2*mt*muf)/B + (16.*mchi2*muf2)/B;
        // g_list[1] = -4.*mchi2 + (8.*mchi*mt)/sqrt(B) + (8.*mchi*muf)/sqrt(B);
        // g_list[2] = 2.;
    }
    else if (oper == 10){
        // g_list[0] = 8*(4*mchi4 - 8*mchi2*mt2 + 4*mt4 - 8*mchi2*s0 - 8*mt2*s0 + 4*s02); 
        // g_list[1] = 8*(-2*mchi2 - 2*mt2 + 4*s0); 
        // g_list[2] = 8.;

        g_list[0] = (-128.*(-1. + B)*mchi2*mt2)/B + (256.*mchi2*mt*muf)/B + (128.*mchi2*muf2)/B;
        g_list[1] = 16.*mchi2 + (64.*mchi*mt)/sqrt(B) + (64.*mchi*muf)/sqrt(B);
        g_list[2] = 8.;
    }
    else if (oper == 11){
        g_list[0] = 0; 
        g_list[1] = 0; 
        g_list[2] = 1;
    }

    double Ci[11] = {2 * mt * mt / vev / vev, 2 * mt * mt / vev / vev, 2 * mt * mt / vev / vev, 2 * mt * mt / vev / vev, 1, 1, 1, 1, 1, 1, 1}; // Operator couplings

    double prefac = sqrt(gamma2) / (2. * s0 * beta - gamma2);
    // double prefac = beta/ (2. * s0 * beta - gamma2)/(16.*pi);
    // double prefac = 0.5*sqrt(1-B)/(mt2+mchi2);

    // double t_min = -(4.*(1-B))/(B*mt2);
    double t_min = -gamma2/s0;
    // double t_min2 = -4.*(1-B)*mt2/B;

    // printf("%.5e\t%0.5e\n", t_min, t_min2);
    // double sum = 0;
    // for (int n=0; n<=2; n++){
    //     sum += g_list[n]*pow(B_fac, n)/(n+1);
    // }
    double sum = 0;
    for (int n = 0; n <= 2; n++)
    {
        // sum += g_list[n] * pow(-1., n)*pow(mt*mchi, 2*(1+n))*pow(mt2 + mchi2, -1-n)*pow(4., 1+n)*(1+mchi2/mt2)*pi*pi/(2*mt2*mchi*(1+n))*(r*r *ne*2.*pow((1-B)/B, n+1)); //* pow(t_min, n) / (n + 1);
        sum += g_list[n] * pow(t_min, n) / (n + 1);// /(n+2);
        // printf("%.5e\t%.5e\", g_list[n], pow(t_min, n));
        // printf("g%d = %0.5e or oper %d with %0.5e\n",n, g_list[n], oper, B);
    }
    // printf("%.5e\n",coeff_high(1));
    

    return prefac*sum*Ci[oper-1];
}

double high_mass_integrand(double x, void * p){

    struct crate_params *params = (struct crate_params *)p;

    double mdm = (params->mchi);
    int np = (params->npts);
    int oper = (params->oper);

    double muf = muFeinterp(x, np);
    double ne = nEinterp(x, np);
    double B = Binterp(x, np);

    // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\n",nstar(B, muf, mdm, oper), mdm, B, muf);
    return x*x*ne*sqrt(1-B)*int_rate_high(mdm, B, muf, oper)/B;///nstar(B, muf, mdm, oper);
    // return int_rate_high(mdm, x, ne, B, muf, oper);
}

double crateHighMass(double mchi, int oper, int npts){

    double res, err;
    size_t cquad_evals;

    struct crate_params params = {mchi, oper, npts};

    size_t sub_divs = 1000;
    // size_t nevals = 20;

    // gsl_integration_workspace *w = gsl_integration_workspace_alloc(nevals);
    // gsl_integration_romberg_workspace *w = gsl_integration_romberg_alloc(nevals);
    gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(sub_divs);

    gsl_function cap_high;

    cap_high.function = &high_mass_integrand;
    cap_high.params = &params;

    // gsl_integration_qags(&cap_high, xmin, 1, 1e-7, 1e-7, nevals,  w, &res, &err);
    // gsl_integration_romberg(&cap_high, xmin, 1, 1e-5, 1e-5, &res, &nevals, w);
    // gsl_integration_qng(&cap_high, xmin, 1, 1e-5, 1e-5, &res, &err, &nevals);
    gsl_integration_cquad(&cap_high, xmin+1e-5, 1, 1e-7, 1e-7, w, &res, &err, &cquad_evals);

    // gsl_integration_workspace_free(w);
    // gsl_integration_romberg_free(w);
    gsl_integration_cquad_workspace_free(w);

    return res*coeff_high(mchi)/2.;


}

//=====================================================
// CAPTURE INTEGRANDS S APPROXIMATION
//=====================================================

double integrandApprx(double r, double uchi, double Ep, double tp, double muFe, double nE, double B, double mchi, int oper, double k0)
{

    double E1 = Eumin(muFe, mchi, k0);
    double E2 = Eumax(muFe);
    double Eu = (E2 - E1) * Ep + E1;

    double delta_S = 4. * sqrt(-1. + (1. + uchi * uchi / 2) * (1. + uchi * uchi / 2) / B) * sqrt(-1. + Eu * Eu) * mt * mchi;
    double s = mt * mt + mchi * mchi + (2. * Eu * mt * mchi * (1. + (uchi * uchi) / 2.)) / sqrt(B);

    double t1 = tmin(s, mchi);
    double t2 = tmax();
    double t = -t1 * tp + t1;

    double EstJac = (E2 - E1) * (t2 - t1) * delta_S; // Jacobian to change (s,t) to unit cube
    // printf("%0.5e\n",tp);
    double soln = solnFunc(Eu, s, t, uchi, mchi, B, muFe);

    double int_main = (r * r * zeta(nE, muFe) * fMB(uchi) * Eu * s * opersdXS(s, t, mchi, oper) * HeavisidePB(Eu, s, t, uchi, mchi, B, muFe, soln) / beta_cap(s, mchi) / gamma_cap(s, mchi) / sqrt(B));

    // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\n", Eu, uchi, B, s1);

    // printf("%f\n", r);
    return EstJac * int_main;
}

//=====================================================
// CAPTURE INTEGRANDS U INTEGRAL APPROXIMATION
//=====================================================

double integrandUappxNoPB(double r, double Ep, double sp, double tp, double muFe, double nE, double B, double mchi, int oper)
{

    // No Pauli Blocking Factor

    double E1 = 1.;
    double E2 = Eumax(muFe);
    double Eu = (E2 - E1) * Ep + E1;

    double s1 = smin(Eu, 0, B, mchi); // smin
    double s2 = smax(Eu, 0, B, mchi); // smax
    double s = (s2 - s1) * sp + s1;

    double t1 = tmin(s, mchi);
    double t2 = tmax();
    double t = -t1 * tp + t1;

    double EstJac = (E2 - E1) * (s2 - s1) * (t2 - t1); // Jacobian to change (Eu, s,t) to unit cube
    // printf("%0.5e\n",tp);
    

    double int_main = (r * r * zeta(nE, muFe) * Eu * s * opersdXS(s, t, mchi, oper) / beta_cap(s, mchi) / gamma_cap(s, mchi) / sqrt(B));
    // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\n", Eu, B, s1);

    return EstJac * int_main;
}

double integrandUappx(double r, double Ep, double sp, double tp, double muFe, double nE, double B, double mchi, int oper, double k0, double z0)
{
    // Full 5D integration

    double E1 = Eumin(muFe, mchi, k0);
    double E2 = Eumax(muFe);
    double Eu = (E2 - E1) * Ep + E1;

    double s1 = smin(Eu, 0, B, mchi); // smin
    double s2 = smax(Eu, 0, B, mchi); // smax
    double s = (s2 - s1) * sp + s1;

    double t1 = (1. - fmax(1. - z0 * mchi / mt, 0))*tmin(s, mchi);
    double t2 = tmax();
    double t = -t1 * tp + t1;

    double EstJac = (E2 - E1) * (s2 - s1) * (t2 - t1); // Jacobian to change (Eu, s, t) to unit cube
    // printf("%0.5e\n",tp);
    double soln = solnFunc(Eu, s, t, 0, mchi, B, muFe);
    // * HeavisidePhaseSpace(Eu, s, B, s1, s2)
    double int_main = (r * r * zeta(nE, muFe) * Eu * s * opersdXS(s, t, mchi, oper) * HeavisidePB(Eu, s, t, 0, mchi, B, muFe, soln)  / beta_cap(s, mchi) / gamma_cap(s, mchi) / sqrt(B));

    // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\n", Eu, B, s1);

    // printf("%f\n", r);
    return EstJac * int_main;
}

double test_integrand(double r, double Eu, double s, double t, double muFe, double nE, double B, double mchi, int oper, double * solout){

    double soln = solnFunc(Eu, s, t, 0, mchi, B, muFe);

    *solout = soln;
    

    double int_main = (r * r * zeta(nE, muFe) * Eu * s * opersdXS(s, t, mchi, oper) * HeavisidePB(Eu, s, t, 0, mchi, B, muFe, soln) / beta_cap(s, mchi) / gamma_cap(s, mchi) / sqrt(B));

    return int_main;
}

//=============================================================================
// MONTE CAPTURE INTEGRANDS U APPRX
//=============================================================================

double monteIntegrandUappxNoPB(double x[], size_t dim, void *p)
{
    (void)(dim);

    struct crate_params *params = (struct crate_params *)p;

    double mdm = (params->mchi);
    int np = (params->npts);
    int oper = (params->oper);


    double r = x[0]; // r/Rstar
    double E = x[1]; // Ep = (Eumax - Eumin)*Eu + Eumin
    double s = x[2]; // sp = (smax - smin)*s + smin
    double t = x[3]; // tp = -tmin*t + tmin

    double muf = muFeinterp(r, np);
    double ne = nEinterp(r, np);
    double B = Binterp(r, np);

    double monteInt = integrandUappxNoPB(r, E, s, t, muf, ne, B, mdm, oper);

    return monteInt;
}

double monteIntegrandUappx(double x[], size_t dim, void *p)
{
    (void)(dim);

    struct crate_params *params = (struct crate_params *)p;

    double mdm = (params->mchi);
    int np = (params->npts);
    int oper = (params->oper);
    double z0 = (params->z0);
    double k0 = (params->k0);

    double r = x[0]; // r/Rstar
    double E = x[1]; // Ep = (Eumax - Eumin)*Eu + Eumin
    double s = x[2]; // sp = (smax - smin)*s + smin
    double t = x[3]; // tp = -tmin*t + tmin

    double muf = muFeinterp(r, np);
    double ne  = nEinterp(r, np);
    double B   = Binterp(r, np);

    double monteInt = integrandUappx(r, E, s, t, muf, ne, B, mdm, oper, k0, z0);

    return monteInt;
}

//=====================================================
// INTEGRATORS U APPROX
//=====================================================

double crateUappxNoPb(double mchi, int oper, int npts)
{

    double res, err;

    struct crate_params params = {mchi, oper, npts};

    double term_max[4] = {1, 1, 1, 1};
    double term_min[4] = {0, 0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function NoPbUappx = {&monteIntegrandUappxNoPB, 4, &params};

    size_t calls = 10000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    {
        gsl_monte_plain_state *s = gsl_monte_plain_alloc(4);
        gsl_monte_plain_integrate(&NoPbUappx, term_min, term_max, 4, calls, rng, s, &res, &err);
        gsl_monte_plain_free(s);

        // display_results("plain", res, err);
    }

    gsl_rng_free(rng);

    return res * CoeffUappx(mchi);
}

double crateUappxFull(double mchi, int oper, int npts, void *cont_vars)
{


    double res, err;

    struct cont_pars *contp  = (struct cont_pars *)cont_vars;

    double z0 = (contp->z0);
    double k0 = (contp->k0);

    struct crate_params params = {mchi, oper, npts, z0, k0};

    double term_max[4] = {1, 1, 1, 1};
    double term_min[4] = {0, 0, 0, 0};

    const gsl_rng_type *Trng;
    gsl_rng *rng;

    gsl_monte_function FullUappx = {&monteIntegrandUappx, 4, &params};

    size_t calls = 50000;

    Trng = gsl_rng_default;
    rng = gsl_rng_alloc(Trng);

    // {
    //     gsl_monte_plain_state *s = gsl_monte_plain_alloc(5);
    //     gsl_monte_plain_integrate(&FullUappx, term_min, term_max, 5, calls, rng, s, &res, &err);
    //     gsl_monte_plain_free(s);

    //     // display_results("plain", res, err);
    // }

    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(4);
        gsl_monte_miser_integrate(&FullUappx, term_min, term_max, 4, calls, rng, s, &res, &err);
        gsl_monte_miser_free(s);

        // display_results("misner", res, err);
        // printf("MISER RESULT:\t%0.5e +- %0.3e\n", res*CoeffUappx(mchi), err*CoeffUappx(mchi));
    }

    // {
    //     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(4);

    //     gsl_monte_vegas_integrate(&FullUappx, term_min, term_max, 4, calls / 5, rng, s, &res, &err);
    //     // display_results("vegas warm-up", res, err);
    //     // printf("converging...\n");

    //     do
    //     {
    //         gsl_monte_vegas_integrate(&FullUappx, term_min, term_max, 4, calls, rng, s, &res, &err);
    //         // printf("result = % .6e sigma = % .6e "
    //         //        "chisq/dof = %.1e\n",
    //         //        res, err, gsl_monte_vegas_chisq(s));
    //     }

    //     while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);

    //     // display_results("vegas final", res, err);

    //     gsl_monte_vegas_free(s);
    //     printf("VEGAS RESULT:\t%0.5e +- %0.3e\n", res, err);
    // }

    gsl_rng_free(rng);

    return res * CoeffUappx(mchi);
}
