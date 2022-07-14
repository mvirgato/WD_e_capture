#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_dawson.h>

#include "globals.h"
#include "collective_effects.h"

/* In this module:
        muFe = Fermi kinetic energy, what we call chemical potential in our early papers
        mue  = Conventional chemical potential, mue = muFe + mt
*/

// Helper functions
double
logabs(double x)
{

    return log(fabs(x));
}

double arcTanh(double x)
{
    return 0.5 * log(fabs((1.0 + x) / (1.0 - x)));
}

double q0_tr(double uchi, double mchi, double B, double muFe, double soln)
{
    // Energy transfer

    double E1 = mchi * (1 + uchi * uchi / 2.) / sqrt(B); // Initial DM energy

    double _q0 = E1 - soln - 0.5 * mchi * uchi * uchi / sqrt(B);

    return _q0;
}

double q_tr_2(double q0, double t)
{

    // Squared 3-momentum transfer

    double _q2 = q0 * q0 - t; // squared momentum transfer

    return _q2;
}

double Eplus(double t, double q0, double q, double mue)
{
    // Kinematic factor (B.8) in Coll. Effects paper

    return fmax(mue, 0.5 * (q0 + q * sqrt(1.0 - 4.0 * mt * mt / t)));
}

// Self Energies for degenerate Fermi gas

double ImPiL_degen(double t, double q0, double q, double muFe)
{
    // Imaginary part of photon self energy

    double mue = muFe + mt; // Change between conventions for muFe. Lasenby uses standard convention compared to our "Fermi Kinetic Energy"
    double ep = Eplus(t, q0, q, mue);

    if (ep > mue + q0)
    {
        return 0.0;
    }

    else
    {
        double ep2 = ep * ep;
        double ep3 = ep * ep * ep;
        double mue2 = mue * mue;
        double mue3 = mue2 * mue;

        double prefac = -echarge * echarge / 4.0 / pi;

        double brak = (-2.0 * ep3 / 3.0 + 2.0 * mue3 / 3.0 + ep2 * q0 + mue2 * q0 - q0 * q0 * q0 / 3.0 - ep * t / 2.0 + mue * t / 2.0 + q0 * t / 2.0);

        return prefac * brak;
    }
}

double ImPiL_degen_approx(double t, double q0, double q, double muFe){

        
    double mue  = muFe + mt; // Change between conventions for muFe. Lasenby uses standard convention compared to our "Fermi Kinetic Energy"
    double pf   = sqrt(mue * mue - mt * mt);
    double Ef   = mue; //sqrt(pf * pf + mt * mt);
    double vf   = pf/Ef;
    double zf = q0 / q / vf;

    double alphaEM = SQR(echarge) / (4.0 * pi);
    double kTF2 = (4.0 * alphaEM / pi) * Ef * pf;

    if (q0 < 0 || fabs(q0) > fabs(q * vf)){
        return 0.0;
    }
    else{

    // return -kTF2 * pi * q0 / 2.0 / vf/ q;
    return - SQR(echarge) * Ef * pf * zf* t /2.0 /pi / SQR(q);
    }
}

double RePiL_degen(double t, double q0, double q, double muFe)
{

    // Real part of photon self energy

    double mue = muFe + mt; // Change between conventions for muFe. Lasenby uses standard convention compared to our "Fermi Kinetic Energy"
    double pf = sqrt(mue * mue - mt * mt);

    double q3 = q * q * q;

    double ik = q * (-(mue * pf) + mt * mt * log((mue + pf) / mt));

    double i0m = integ0m(t, q0, q, mue, pf);
    double i1m = integ1m(t, q0, q, mue, pf);
    double i2m = integ2m(t, q0, q, mue, pf);
    double i0p = integ0p(t, q0, q, mue, pf);
    double i1p = integ1p(t, q0, q, mue, pf);
    double i2p = integ2p(t, q0, q, mue, pf);

    double res = (ik - (t / 4) * i0m - i2m + q0 * i1m + i2p + (t / 4) * i0p + q0 * i1p);
    return (t / q3) * (echarge * echarge / (2 * pi * pi)) * res;
}

double RePiL_degen_approx(double t, double q0, double q, double muFe)
{

    // Approximate Real part of photon self energy for small momentum transfer

    double mue = muFe + mt; // Change between conventions for muFe. Lasenby uses standard convention compared to our "Fermi Kinetic Energy"
    double pf = sqrt(mue * mue - mt * mt);
    double Ef = mue; //sqrt(pf * pf + mt * mt);
    double vf = pf / Ef;
    double zf = q0 / q / vf;

    double prefac = SQR(echarge) / SQR(pi);

    double res = prefac * (t / SQR(q)) * Ef * pf * (-1.0 + 0.5 * zf * logabs((1 + zf) / (1 - zf)));

    return res;
}

double RePiL_degen_approx_2(double t, double q0, double q, double muFe){

    double mue  = muFe + mt; // Change between conventions for muFe. Lasenby uses standard convention compared to our "Fermi Kinetic Energy"
    double pf   = sqrt(mue * mue - mt * mt);
    double Ef   = mue; //sqrt(pf * pf + mt * mt);
    // double vf   = pf/Ef;

    double alphaEM = SQR(echarge) / (4.0 * pi);
    double kTF2 = (4.0 *  alphaEM / pi) * Ef * pf;

    return kTF2;
}

// double PiL_degen(double t, double uchi, double mchi, double B, double muFe, double soln){

//     // Full longintudinal self energy

//     double q0 = q0_tr(uchi, mchi, B, muFe, soln);
//     double q  = sqrt(q_tr_2(q0, t));

//     return

// }


// Dilute non-relativistic gas

double Zi(double x)
{

    return sqrt(pi) * exp(-x * x);
}

double Zr(double x)
{

    return -sqrt(pi) * gsl_sf_dawson(x);
}

double RePiL_dnr(double t, double q0, double q, double nE)
{
    double m = 12000.0;
    nE = nE / (pmTOm * mTOinveV)/ (pmTOm * mTOinveV)/ (pmTOm * mTOinveV);
    double T = 1e5 * kBMeV;
    double Z = 6.0;
    double s = sqrt(T / m);
    double xi = q0 / (sqrt(2) * q * s);
    double qp = t / (2 * q0 * m);
    return (((SQR(echarge) * SQR(Z) * nE/Z) / (sqrt(2) * q * s)) * (Zr(xi * (1.0 + qp)) - Zr(xi * (1.0 - qp))));
}

double ImPiL_dnr(double t, double q0, double q, double nE)
{
    double m = 12000.0;
    nE = nE / (pmTOm * mTOinveV)/ (pmTOm * mTOinveV)/ (pmTOm * mTOinveV);
    double T = 1e5 * kBMeV;
    double Z = 6.0;
    double s = sqrt(T / m);
    double xi = q0 / (sqrt(2) * q * s);
    double qp = t / (2 * q0 * m);

    return (((SQR(echarge) * SQR(Z) * nE/Z) / (sqrt(2) * q * s)) * (Zi(xi * (1.0 + qp)) - Zi(xi * (1.0 - qp))));
}

double RePiL_dnr_approx(double t, double q0, double q, double nE){

    double mi  = 12000.0;
    nE = nE / (pmTOm * mTOinveV)/ (pmTOm * mTOinveV)/ (pmTOm * mTOinveV);
    double wi2 = (nE/6.0) * SQR(6.0) * SQR(echarge) / mi;
    

    return wi2 * t / SQR(q0);
}



double PiL_appox(double t, double q0, double q, double muFe, double nE){

    double mi  = 12000.0;
    nE = nE / (pmTOm * mTOinveV)/ (pmTOm * mTOinveV)/ (pmTOm * mTOinveV);
    double wi2 = (nE/6.0) * SQR(6.0) * SQR(echarge) / mi;

    double mue = muFe + mt; // Change between conventions for muFe. Lasenby uses standard convention compared to our "Fermi Kinetic Energy"
    double pf = sqrt(mue * mue - mt * mt);
    double Ef = mue; //sqrt(pf * pf + mt * mt);

    double alphaEM = SQR(echarge) / (4.0 * pi);

    double PiIon  = wi2 * t / SQR(q0);
    double PiELec = - (4.0 *  alphaEM / pi) * Ef * pf * t/SQR(q);

    return PiELec + PiIon;
}

double CollEffectsFF(double t, double uchi, double mchi, double B, double muFe, double nE, double soln)
{
    if (t == 0){
        return 0.0;
    }
    else{

        double q0 = q0_tr(uchi, mchi, B, muFe, soln);
        double q = sqrt(q_tr_2(q0, t));

        // Fudge Factor to account for collective effects

        double RePiL = RePiL_degen_approx(t, q0, q, muFe) + RePiL_dnr(t, q0, q, nE);
        double ImPiL = ImPiL_degen_approx(t, q0, q, muFe) + ImPiL_dnr(t, q0, q, nE);

        // printf("%0.5e\n", RePiL_degen_approx(t, q0, q, muFe) / RePiL_degen(t, q0, q, muFe));

        double denom = ((RePiL - t) * (RePiL - t) + ImPiL * ImPiL);
        // double denom = SQR(t - PiL_appox(t, q0, q, muFe, nE)) + PiL_appox(t, q0, q, muFe, nE)/10.0;

        double res = t * t / denom;
        // if (res > 1.0){
        //     printf("%0.5e\n", res);
        // }

        return res;
    }
}

// Analytic expressions for integrals

double integ0m(double t, double q0, double q, double muFe, double pf)
{

    double QQ2 = q0 * q0 + q * q;

    double C1 = sqrt(t * (t - 4.0 * mt * mt));

    return ((2 * q * logabs((muFe + pf) / mt) + 2 * muFe * logabs((-2 * pf * q - 2 * muFe * q0 + t) / (2 * pf * q - 2 * muFe * q0 + t)) + ((4 * (mt * mt) * q * q0 - 2 * q * q0 * t + C1 * QQ2) * (logabs(mt) + logabs(-(C1 * q) - 2 * muFe * t + q0 * t))) / sqrt(t * (-2 * q0 * (C1 * q + 2 * (mt * mt) * q0) + t * QQ2)) + ((2 * q * q0 * (-2 * (mt * mt) + t) + C1 * QQ2) * (-logabs(mt) - logabs(C1 * q - 2 * muFe * t + q0 * t))) / sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)) - ((4 * (mt * mt) * q * q0 - 2 * q * q0 * t + C1 * QQ2) * logabs(C1 * muFe * q + 2 * (mt * mt) * t - muFe * q0 * t + pf * sqrt(t * (-2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)))) / sqrt(t * (-2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)) + ((-4 * (mt * mt) * q * q0 + 2 * q * q0 * t + C1 * QQ2) * logabs(C1 * muFe * q - 2 * (mt * mt) * t + muFe * q0 * t - pf * sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)))) / sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2))) / 2);
}

double integ0p(double t, double q0, double q, double muFe, double pf)
{

    double QQ2 = q0 * q0 + q * q;

    double C1 = sqrt(t * (t - 4 * (mt * mt)));

    return ((-2 * q * logabs((muFe + pf) / mt) + 2 * muFe * logabs((2 * pf * q + 2 * muFe * q0 + t) / (-2 * pf * q + 2 * muFe * q0 + t)) + ((4 * (mt * mt) * q * q0 - 2 * q * q0 * t - C1 * QQ2) * (logabs(mt) + logabs(-(C1 * q) - (2 * muFe + q0) * t))) / sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)) - ((4 * (mt * mt) * q * q0 - 2 * q * q0 * t + C1 * QQ2) * (-logabs(mt) - logabs(C1 * q - (2 * muFe + q0) * t))) / sqrt(t * (-2 * q0 * (C1 * q + 2 * (mt * mt) * q0) + t * QQ2)) - ((4 * (mt * mt) * q * q0 - 2 * q * q0 * t + C1 * QQ2) * logabs(C1 * muFe * q - 2 * (mt * mt) * t - muFe * q0 * t - pf * sqrt(t * (-2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)))) / sqrt(t * (-2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)) + ((-4 * (mt * mt) * q * q0 + 2 * q * q0 * t + C1 * QQ2) * logabs(C1 * muFe * q + 2 * (mt * mt) * t + muFe * q0 * t + pf * sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)))) / sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2))) / 2);
}

double integ1m(double t, double q0, double q, double muFe, double pf)
{

    double QQ2 = q0 * q0 + q * q;

    double C1 = sqrt(t * (t - 4 * (mt * mt)));

    return ((pf * q + (q * q0 * (-2 * (mt * mt) + t) * arcTanh(muFe / pf)) / t - ((8 * q0 * t * (-2 * (mt * mt) + t) * (4 * (mt * mt) * (q * q) + (t * t)) + (4 * C1 * q - 4 * q0 * t) * (4 * (mt * mt) * ((q * q * q * q) + (q * q) * (q0 * q0) - 2 * (q0 * q0 * q0 * q0)) + ((q * q) + 3 * (q0 * q0)) * (t * t))) * sqrt(t * (-2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)) * arcTanh((-(C1 * muFe * q) - 2 * (mt * mt) * t + muFe * q0 * t) / (sqrt(-(mt * mt) + (muFe * muFe)) * sqrt(t * (-2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2))))) / (C1 * t * (-64 * (mt * mt) * (t * t) + (4 * C1 * q - 4 * q0 * t) * (4 * C1 * q - 4 * q0 * t))) - ((16 * (mt * mt * mt * mt) * q0 * ((q * q * q) - q * (q0 * q0)) * ((q * q * q) - q * (q0 * q0)) + (t * t * t) * (C1 * ((q * q * q) + 3 * q * (q0 * q0)) + q0 * (3 * (q * q) + (q0 * q0)) * t) - 4 * (mt * mt) * (t * t) * (C1 * ((q * q * q) + 2 * q * (q0 * q0)) + q0 * (4 * (q * q) + (q0 * q0)) * t)) * arcTanh((C1 * muFe * q - 2 * (mt * mt) * t + muFe * q0 * t) / (sqrt(-(mt * mt) + (muFe * muFe)) * sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2))))) / (4 * C1 * (t * t) * sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2))) + (muFe * muFe) * logabs((-2 * pf * q - 2 * muFe * q0 + t) / (2 * pf * q - 2 * muFe * q0 + t))) / 2);
}

double integ1p(double t, double q0, double q, double muFe, double pf)
{

    double QQ2 = q0 * q0 + q * q;

    double C1 = sqrt(t * (t - 4 * (mt * mt)));

    return ((-(pf * q) + (q * q0 * (-2 * (mt * mt) + t) * arcTanh(muFe / pf)) / t + ((8 * q0 * t * (-2 * (mt * mt) + t) * (4 * (mt * mt) * (q * q) + (t * t)) - (-4 * C1 * q + 4 * q0 * t) * (4 * (mt * mt) * ((q * q * q * q) + (q * q) * (q0 * q0) - 2 * (q0 * q0 * q0 * q0)) + ((q * q) + 3 * (q0 * q0)) * (t * t))) * sqrt(t * (-2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)) * arcTanh((C1 * muFe * q - (2 * (mt * mt) + muFe * q0) * t) / (pf * sqrt(t * (-2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2))))) / (C1 * t * (-64 * (mt * mt) * (t * t) + (4 * C1 * q - 4 * q0 * t) * (4 * C1 * q - 4 * q0 * t))) + ((16 * (mt * mt * mt * mt) * q0 * ((q * q * q) - q * (q0 * q0)) * ((q * q * q) - q * (q0 * q0)) + t * t * t * (C1 * ((q * q * q) + 3 * q * (q0 * q0)) + q0 * (3 * (q * q) + (q0 * q0)) * t) - 4 * (mt * mt) * (t * t) * (C1 * ((q * q * q) + 2 * q * (q0 * q0)) + q0 * (4 * (q * q) + (q0 * q0)) * t)) * arcTanh((-2 * (mt * mt) * t - muFe * (C1 * q + q0 * t)) / (pf * sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2))))) / (4 * C1 * (t * t) * sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2))) + (muFe * muFe) * logabs((2 * pf * q + 2 * muFe * q0 + t) / (-2 * pf * q + 2 * muFe * q0 + t))) / 2);
}

double integ2m(double t, double q0, double q, double muFe, double pf)
{

    double QQ2 = q0 * q0 + q * q;

    double C1 = sqrt(t * (t - 4 * (mt * mt)));

    return ((2 * muFe * pf * q + (4 * pf * q * q0 * (-2 * (mt * mt) + t)) / t + 2 * (mt * mt) * q * arcTanh(muFe / pf) + (q * (4 * (mt * mt) * ((q * q * q * q) + (q * q) * (q0 * q0) - 2 * (q0 * q0 * q0 * q0)) + ((q * q) + 3 * (q0 * q0)) * (t * t)) * arcTanh(muFe / pf)) / (t * t) - (8 * sqrt(t * (-2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)) * ((t * t) * (4 * (mt * mt) * (q * q) + (t * t)) * (4 * (mt * mt) * ((q * q * q * q) + (q * q) * (q0 * q0) - 2 * (q0 * q0 * q0 * q0)) + ((q * q) + 3 * (q0 * q0)) * (t * t)) + q0 * (4 * C1 * q - 4 * q0 * t) * (4 * (mt * mt * mt * mt) * ((q * q * q) - q * (q0 * q0)) * ((q * q * q) - q * (q0 * q0)) + (mt * mt) * (5 * (q * q * q * q) - 2 * (q * q) * (q0 * q0) - 3 * (q0 * q0 * q0 * q0)) * (t * t) + (t * t * t * t) * QQ2)) * arcTanh((-(C1 * muFe * q) - 2 * (mt * mt) * t + muFe * q0 * t) / (sqrt(-(mt * mt) + (muFe * muFe)) * sqrt(t * (-2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2))))) / (C1 * t * t * t * (-64 * (mt * mt) * (t * t) + (4 * C1 * q - 4 * q0 * t) * (4 * C1 * q - 4 * q0 * t))) + (sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)) * ((t * t) * (4 * (mt * mt) * (q * q) + (t * t)) * (4 * (mt * mt) * ((q * q * q * q) + (q * q) * (q0 * q0) - 2 * (q0 * q0 * q0 * q0)) + ((q * q) + 3 * (q0 * q0)) * (t * t)) - 4 * q0 * (C1 * q + q0 * t) * (4 * (mt * mt * mt * mt) * ((q * q * q) - q * (q0 * q0)) * ((q * q * q) - q * (q0 * q0)) + (mt * mt) * (5 * (q * q * q * q) - 2 * (q * q) * (q0 * q0) - 3 * (q0 * q0 * q0 * q0)) * (t * t) + (t * t * t * t) * QQ2)) * arcTanh((C1 * muFe * q - 2 * (mt * mt) * t + muFe * q0 * t) / (sqrt(-(mt * mt) + (muFe * muFe)) * sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2))))) / (2 * C1 * t * t * t * (-4 * (mt * mt) * (t * t) + (C1 * q + q0 * t) * (C1 * q + q0 * t))) + 4 * muFe * muFe * muFe * logabs((-2 * pf * q - 2 * muFe * q0 + t) / (2 * pf * q - 2 * muFe * q0 + t))) / 12);
}

double integ2p(double t, double q0, double q, double muFe, double pf)
{

    double QQ2 = q0 * q0 + q * q;

    double C1 = sqrt(t * (t - 4 * (mt * mt)));
    return ((-2 * muFe * pf * q + (4 * pf * q * q0 * (-2 * (mt * mt) + t)) / t - 2 * (mt * mt) * q * arcTanh(muFe / pf) - (q * (4 * (mt * mt) * ((q * q * q * q) + (q * q) * (q0 * q0) - 2 * (q0 * q0 * q0 * q0)) + ((q * q) + 3 * (q0 * q0)) * (t * t)) * arcTanh(muFe / pf)) / (t * t) - (sqrt(t * (-2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)) * (8 * (t * t) * (4 * (mt * mt) * (q * q) + (t * t)) * (4 * (mt * mt) * ((q * q * q * q) + (q * q) * (q0 * q0) - 2 * (q0 * q0 * q0 * q0)) + ((q * q) + 3 * (q0 * q0)) * (t * t)) - 8 * q0 * (-4 * C1 * q + 4 * q0 * t) * (4 * (mt * mt * mt * mt) * ((q * q * q) - q * (q0 * q0)) * ((q * q * q) - q * (q0 * q0)) + (mt * mt) * (5 * (q * q * q * q) - 2 * (q * q) * (q0 * q0) - 3 * (q0 * q0 * q0 * q0)) * (t * t) + (t * t * t * t) * QQ2)) * arcTanh((C1 * muFe * q - (2 * (mt * mt) + muFe * q0) * t) / (pf * sqrt(t * (-2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2))))) / (C1 * t * t * t * (-64 * (mt * mt) * (t * t) + (4 * C1 * q - 4 * q0 * t) * (4 * C1 * q - 4 * q0 * t))) + (sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2)) * (8 * (t * t) * (4 * (mt * mt) * (q * q) + (t * t)) * (4 * (mt * mt) * ((q * q * q * q) + (q * q) * (q0 * q0) - 2 * (q0 * q0 * q0 * q0)) + ((q * q) + 3 * (q0 * q0)) * (t * t)) - 32 * q0 * (C1 * q + q0 * t) * (4 * (mt * mt * mt * mt) * ((q * q * q) - q * (q0 * q0)) * ((q * q * q) - q * (q0 * q0)) + (mt * mt) * (5 * (q * q * q * q) - 2 * (q * q) * (q0 * q0) - 3 * (q0 * q0 * q0 * q0)) * (t * t) + (t * t * t * t) * QQ2)) * arcTanh((-2 * (mt * mt) * t - muFe * (C1 * q + q0 * t)) / (pf * sqrt(t * (2 * C1 * q * q0 - 4 * (mt * mt) * (q0 * q0) + t * QQ2))))) / (16 * C1 * t * t * t * (-4 * (mt * mt) * (t * t) + (C1 * q + q0 * t) * (C1 * q + q0 * t))) + 4 * muFe * muFe * muFe * logabs((2 * pf * q + 2 * muFe * q0 + t) / (-2 * pf * q + 2 * muFe * q0 + t))) / 12);
}