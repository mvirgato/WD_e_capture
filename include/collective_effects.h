#ifndef COLLEFF_H
#define COLLEFF_H

double integ0m(double t, double q0, double q, double muFe, double pf);
double integ1m(double t, double q0, double q, double muFe, double pf);
double integ2m(double t, double q0, double q, double muFe, double pf);
double integ0p(double t, double q0, double q, double muFe, double pf);
double integ1p(double t, double q0, double q, double muFe, double pf);
double integ2p(double t, double q0, double q, double muFe, double pf);

double ImPiL_degen(double t, double q0, double q, double muFe);
double RePiL_degen(double t, double q0, double q, double muFe);
double CollEffectsFF(double t, double uchi, double mchi, double B, double muFe, double nE, double soln);

double Zr(double x);
double RePiL_dnr(double t, double q0, double q, double nE);
double ImPiL_dnr(double t, double q0, double q, double nE);

#endif