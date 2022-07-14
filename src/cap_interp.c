#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "globals.h"
#include "cap_interp.h"



#define numprofs 5000

double xint[numprofs];
double Mint[numprofs];
double densint[numprofs];
double Pint[numprofs];
double nEint[numprofs];
double nPint[numprofs];
double muFeint[numprofs];
double Bint[numprofs];

double Bms[100];
double mufms[100];
double mstar[10][100][100];
int num_B  = 100;
int num_mu = 100;
double step_B = 0.0005049494949494959;
double step_mu = 0.21;


int read_profs(char *filename, int *np)
{
    // reads in EOS radial profiles for B, muFe, and nE as functions of r/Rstar (0, 1)
    // order: nE (pm^-3); muFe (MeV); B 

    FILE *datafile;
    int i, npts;

    datafile = fopen(filename, "r");

    if (datafile == NULL)
    {
        printf("Error: can't open file %s\n", filename);
        exit(0);
        return 1;
    }

    else
    {
        printf("File %s opened successfully.\n", filename);
        i = 0;
        fscanf(datafile, "%*[^\n]\n", NULL);
        while (!feof(datafile))
        {
            if(fscanf(datafile, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &xint[i], &Mint[i], &densint[i], &Pint[i], &nEint[i], &nPint[i], &muFeint[i], &Bint[i]));
            // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n", xint[i], Mint[i], densint[i], Pint[i], nEint[i], muFeint[i], Bint[i]);
            i++;
        }
    }

    fclose(datafile);

    RS = xint[numprofs-1];
    xmin = xint[0]/RS;
    rhoc = densint[0];

    int j;
    for (j = 0; j<=numprofs-1; j++){
        xint[j] = xint[j]/RS;
    }


    *np = i;
    n0 = nEint[0];
    muFe0 = muFeint[0];
    B0 = Bint[0];

    return 0;
}

int read_mstar(char *filename)
{   
    double Btemp[10000];
    double mutemp[10000];
    double mstemp[10][10000];



    FILE *datafile;
    int i, npts;

    datafile = fopen(filename, "r");

    if (datafile == NULL)
    {
        printf("Error: can't open file %s\n", filename);
        exit(0);
        return 1;
    }

    else
    {
        printf("File %s opened successfully.\n", filename);
        i = 0;
        fscanf(datafile, "%*[^\n]\n", NULL);
        while (!feof(datafile))
        {
            if (fscanf(datafile, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                       &Btemp[i], &mutemp[i], &mstemp[0][i], &mstemp[1][i], &mstemp[2][i], &mstemp[3][i], &mstemp[4][i], &mstemp[5][i], &mstemp[6][i], &mstemp[7][i], &mstemp[8][i], &mstemp[9][i]));
                // printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                //        Btemp[i], mutemp[i], mstemp[0][i], mstemp[1][i], mstemp[2][i], mstemp[3][i], mstemp[4][i], mstemp[5][i], mstemp[6][i], mstemp[7][i], mstemp[8][i], mstemp[9][i]);
            i++;
        }
    }
    fclose(datafile);

    
    double bmin = Btemp[0];
    double bmax = Btemp[i-1];

    double mumin = mutemp[0];
    double mumax = mutemp[i-1];

    int l, m;
    for (l = 0; l < 100; l++)
    {
        Bms[l] = bmin + step_B * l;
        // printf("%f\n", Bms[l]);
    }

    for (m = 0; m < 100; m++)
    {
        mufms[m] = mumin + step_mu * m;
        // printf("%f\n", mufms[m]);
    }

    int op;
    for (op = 0; op < 10; op++ ){
        for (l = 0; l < 100; l++){
            for (m = 0; m < 100; m++){
                mstar[op][l][m] = mstemp[op][l*100 + m];
            }
        }
    }

    return 0;
}

double Binterp(double x, int npts)
{

    double B_r;

    if (x >= xint[0] && x <= xint[npts - 1])
    {

        gsl_interp_accel *acc = gsl_interp_accel_alloc();

        // Cubic splines
        gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, npts);

        gsl_spline_init(spline, xint, Bint, npts);

        B_r = gsl_spline_eval(spline, x, acc);

        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);

        return B_r;
    }

    else
    {
        return 0.;
    }
}

double muFeinterp(double x, int npts)
{

    double muFe_r;

    if (x >= xint[0] && x <= xint[npts - 1])
    {

        gsl_interp_accel *acc = gsl_interp_accel_alloc();

        // Cubic splines
        gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, npts);

        gsl_spline_init(spline, xint, muFeint, npts);

        muFe_r = gsl_spline_eval(spline, x, acc);

        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);

        return muFe_r;
    }

    else
    {
        return 0.;
    }
}

double nEinterp(double x, int npts)
{

    double nE_r;

    if (x >= xint[0] && x <= xint[npts - 1])
    {

        gsl_interp_accel *acc = gsl_interp_accel_alloc();

        // Cubic splines
        gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, npts);

        gsl_spline_init(spline, xint, nEint, npts);

        nE_r = gsl_spline_eval(spline, x, acc);

        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);


        return nE_r;
    }

    else
    {
        return 0.;
    }
}

double mstarInterp(double Bin, double mufin, int oper)
{
    int npts = 100;
    if (Bin >= Bms[0] && Bin <= Bms[npts - 1] && mufin >= mufms[0] && mufin <= mufms[npts - 1])
    {

        int ms_num = 100*100;
        double ms_arr[ms_num];

        double ms_r;

        int i, j;
        for (i = 0; i<npts; i++){
            for (j = 0; j < npts; j++){
                ms_arr[j*npts + i] = mstar[oper-1][i][j];
            }
        }

        const gsl_interp2d_type *T = gsl_interp2d_bicubic;

        gsl_spline2d *spline = gsl_spline2d_alloc(T, npts, npts);
        gsl_interp_accel *bacc = gsl_interp_accel_alloc();
        gsl_interp_accel *muacc = gsl_interp_accel_alloc();

        gsl_spline2d_init(spline, Bms, mufms, ms_arr, npts, npts);

        ms_r = gsl_spline2d_eval(spline, Bin, mufin, bacc, muacc);
        // printf("%0.5f\n", interp_lFval);

        gsl_spline2d_free(spline);
        gsl_interp_accel_free(bacc);
        gsl_interp_accel_free(muacc);



        return ms_r;
    }

    else
    {
        printf("otside of range:\n");
        if (Bin < Bms[0] || Bin > Bms[npts-1]){
            printf(" Bin = %f\n", Bin);
        }
        if (mufin < mufms[0] || mufin > mufms[npts-1]){
            printf(" mufin = %f\n", mufin);
        }
            return 1.;
    }
}