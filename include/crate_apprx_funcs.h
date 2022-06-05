#ifndef CAPAPPROX_H
#define CAPAPPROX_H

// General Funcs
double CoeffUappx(double mchi);
double int_rate_high(double mchi, double B, double muf, int oper);
// double int_rate_high(double mchi, double r, double ne, double B, double muf, int oper);
double crateHighMass(double mchi, int oper, int npts);
double cap_high_2(double mchi, int oper, int np);

// Integrands
double ds10r(double mchi, double muf, double B, int n, int oper);
// double integrandApprx(double r, double uchi, double Ep, double tp, double muFe, double nE, double B, double mchi, int oper);
double integrandUappx(double r, double Ep, double sp, double tp, double muFe, double nE, double B, double mchi, int oper, double k0, double z0);
double integrandUappxNoPB(double r, double Ep, double sp, double tp, double muFe, double nE, double B, double mchi, int oper);
double test_integrand(double r, double Eu, double s, double t, double muFe, double nE, double B, double mchi, int oper, double * solout);

// Integrators
double crateUappxFull(double mchi, int oper, int npts, void *cont_vars);
double crateUappxNoPb(double mchi, int oper, int npts);

#endif