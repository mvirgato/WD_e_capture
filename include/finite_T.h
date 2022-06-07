#ifndef finit_T_H
#define finit_T_H

double Coeff_temp(double mchi, double temp);
double CoeffUappx_temp(double mchi, double temp);
double crateNoPb_T(double mchi, double temp, int oper, int npts, void *cont_vars);
double crateFull_T(double mchi, double temp, int oper, int npts, void *cont_vars);
double crateUappx_T(double mchi, double temp, int oper, int npts, void *cont_vars);

#endif