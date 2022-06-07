#ifndef evap_func_H
#define evap_func_H

//misc functions

// evap functions
double radChiSqr(double T, double mchi);
double niso(double rad, double T, double mchi);
double fMBevap(double vel, double vesc, double T, double mchi);
double prefactorEvap(double mchi, double T);
double heavisideprod(double w, double v, double t, double s);

// Integrands
double evapIntegrand(double r, double wp, double vp, double tp, double sp, double mchi, double vesc, double T, double ne);

// Monte Integrands
double monteEvapIntegrand(double x[], size_t dim, void * p);

// Integrators
struct evap_params {double mchi; double T; int npts;};
double evaprate(double mchi, double T, int npts);
double reg1evap(double mchi, double T, int npts);

#endif