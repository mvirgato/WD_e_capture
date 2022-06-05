#ifndef CAPINT_H
#define CAPINT_H

double n0; // central electron number density
double muFe0; // central electron chemical potential
double RS; // Star radius
double xmin;// min x = r/RS
double rhoc; // central density

int read_profs(char *filename, int *np);
int read_mstar(char *filename);
double Binterp(double x, int npts);
double muFeinterp(double x, int npts);
double nEinterp(double x, int npts);
double mstarInterp(double Bin, double mufin, int oper);

#endif