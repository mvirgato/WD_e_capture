#ifndef CAPRAD_H
#define CAPRAD_H

double crateRaProfFull(double rad, double mchi, int oper, int npts);
double crateRaProfNoPb(double rad, double mchi, int oper, int npts);
double crateRaProfAppx(double rad, double mchi, int oper, int npts, int iters);

#endif
