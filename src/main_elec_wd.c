#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>

#include <omp.h>

//=============================================================================

#include "globals.h"
#include "cap_interp.h"
#include "cap_funcs.h"
#include "crate_apprx_funcs.h"
// #include "crate_radial_profs.h"
#include "finite_T.h"
#include "multiscatter.h"
// #include "use_cuba.h"
#include "collective_effects.h"
#include "interaction_rates.h"

//=============================================================================

int np;
char *set_WD;
char *elem;
char *label;

struct cont_pars control_params;

char *WD_names[44] = {"0.304035", "0.459362", "0.49", "0.569974", "0.67", "0.960626", "1", "1.04032", "1.15", "1.18045", "1.19646", "1.20678", "1.21022", "1.21776", "1.22043", "1.22278", "1.23266", "1.23281", "1.24103", "1.24953", "1.24985", "1.25278", "1.25896", "1.26003", "1.28434", "1.28641", "1.29718", "1.30287", "1.31879", "1.32347", "1.32525", "1.33188", "1.33559", "1.34482", "1.35776", "1.36", "1.3637", "1.36814", "1.37195", "1.37215", "1.38067", "1.3818", "1.38222", "1.38367"};

// char *WD_names[44] = {"0.598985","0.628639","0.632971","1.24138","1.31033","1.61194","1.617","1.71238","1.80496","1.92247","2.146","2.29479","2.34608","2.43769","2.46225","2.52696","2.72818","2.79624","2.92496","2.94971","3.20695","3.21817","3.283","3.31376","3.31715","3.40178","3.479","3.48038","3.57319","3.59525","3.62036","3.69112","3.72351","3.81837","3.95526","4.2177","5.07458","5.37714","5.66377","7.81997","8.6581","9.38577","9.69494","11.5929"};
char * elem_list[3] = {"He", "C", "O"};

int USE_SPARTAN; // set = 1 if using spartan (or other clusters), 0 if not

//=============================================================================
// void method_comps();
void scan_1D_region(double mchi, int oper, int start, int end, float step);
void scan_2D_region(double mchi, int oper, int start, int end, float step);

double crate_complete(double lmmin, double lmmax, double step);

double crate_singleOp_Full(int oper, double lmmin, double lmmax, double step);
double crate_singleOp_CollEff(int oper, double lmmin, double lmmax, double step);
double crate_singleOp_NoPB(int oper, double lmmin, double lmmax, double step);
double crate_singleOp_Uappx_Full(int oper, double lmmin, double lmmax, double step);
double crate_all_Full(double lmmin, double lmmax, double step);
double crate_all_NoPB(double lmmin, double lmmax, double step);
double crate_all_Uappx_Full(double lmmin, double lmmax, double step);
double crate_all_Uappx_NoPB(double lmmin, double lmmax, double step);

double int_rate_singleOp_Full(double r, int oper, double lmmin, double lmmax, double step);

double crate_singleOp_finiteT_full(int oper, double lmmin, double lmmax, double temp, double step);
double crate_all_finiteT_NoPB(double lmmin, double lmmax, double temp, double step);
double crate_all_finiteT_Full(double lmmin, double lmmax, double temp, double step);
double crate_all_finiteT_Uappx(double lmmin, double lmmax, double temp, double step);

double crate_all_multiscat(double lmmin, double lmmax, double step);

double cap_prof_full(double rmin, double rmax, double mchi, int num_pts);
double cap_prof_nopb(int oper, double rmin, double rmax, double mchi, double step);
double cap_prof_appx(int oper, double rmin, double rmax, double mchi, double step);

//=============================================================================
//=============================================================================
// BEGIN MAIN
//=============================================================================
//=============================================================================

int main(int argc, char *argv[])
{

//=============================================================================
// INITIALIZATION FROM COMMAND LINE INPUTS AND USER INPUTS
//=============================================================================

elem = elem_list[1]; // set element here: 0 = He, 1 = C, 2 = O
label = "M"; //set to "M" for mass labels, or "r" for radius labels

if (argc == 1)
{
        set_WD = "1.38367"; // Change this to set a different WD
        printf("No WD selected.\nUsing default %s WD with %s EoS\n", set_WD, elem);
        USE_SPARTAN = 0;
}
else{
        printf("Selected %s M_sun WD with %s EoS.\n", WD_names[atoi(argv[1])], elem);
        set_WD = WD_names[atoi(argv[1])];
        USE_SPARTAN = atoi(argv[2]);
}


if (USE_SPARTAN == 0)
{
        int NUM_CORES = omp_get_max_threads(); // Default number of cores when not useing Spartan
        omp_set_num_threads(NUM_CORES);
        printf("Running locally on %d cores\n", NUM_CORES);
}

system_setup(mE, set_WD, label, elem, "20", &np); //outputs mt = target mass and loads profiles

// control parameters

control_params.z0 = 0.005;
control_params.k0 = 0.0005;
control_params.ct = 1.;



//=============================================================================
// CAPTURE CALCULATION
//=============================================================================

        // crate_singleOp_NoPB(5, 0, 4, 0.2);
        // crate_singleOp_Full(5, 0., 4., 0.2);
        // crate_singleOp_Uappx_Full(1, -3, 3., 0.2);
        // crate_singleOp_CollEff(5, 0, 4, 0.2);
        // if (fabs(x))
        //  crate_all_NoPB(2., 8., 0.2);
        //  crate_all_Full(-2., 8., 0.2);

        //  crate_all_Uappx_NoPB(-2., 8., 0.2);
        //  crate_all_Uappx_Full(-2., 8., 0.2);


        // crate_all_multiscat(3., 6., 0.2);

        // crate_complete(-2, 8, 0.2);

        // double soln = solnFunc(0.1, 0.2, -0.001, 1e-3, 0.1, 0.5, 10.0);
        // printf("%0.5e\n", soln);
        // double out = CollEffectsFF(-0.001, 1e-3, 0.1, 0.5, 10.0, soln);

        // printf("%0.5e\n", out);

        //=============================================================================
        // INTERACTION RATES CALCULATION
        //=============================================================================
        // printf("%0.5e\n", xmin);
        int_rate_singleOp_Full(xmin, 5, 0, 3, 0.2);
        //=============================================================================
        // RADIAL PROFILE CALCULATION
        //=============================================================================

        // cap_prof_nopb(11,  0.9, 1, 0.1, 0.005);
        // cap_prof_full(0.001, 1., 1., 50);

        // cap_prof_appx(11, 0., 1, 1, 0.1);

        //=============================================================================
        // FINITE TEMPERATURE CALCULATION
        //=============================================================================

        // double Tstar = kBMeV * 1e5;
        // crate_all_finiteT_Full(-2., 8., Tstar, 0.2);
        // crate_all_finiteT_Uappx(-2., 8, Tstar, 0.2);


//=============================================================================
// TESTS
//=============================================================================
        // printf("%0.5e\n", crateMultiScatt(1e5, 1, np));
        // int oper = 4;
        // int mchi = 1;

        // double out1 = crateUappxFull(mchi, oper, np, &control_params);
        // double out2 = crateUappx_T(mchi, Tstar, oper, np, &control_params);

        // double out1 = crateFull(mchi, oper, np, &control_params);
        // double out2 = crateFull_T(mchi, Tstar, oper, np, &control_params);

        // printf("%f\n", out2out1);

        // double out = Zr(1);
        // printf("%0.5e\n", out);
//=============================================================================
//=============================================================================

return 0;
}

// END MAIN
//=============================================================================
//=============================================================================


//=============================================================================
// INTERACTION RATES
//=============================================================================

double int_rate_singleOp_Full(double r, int oper, double lmmin, double lmmax, double step)
{
        // enter r in r/Rstar

        int i; // mass loop
        int end = round((lmmax - lmmin) / step)+1; //number of mass points
        double intOut[end];
        double mchi[end];
        int tally = 0;

        #pragma omp parallel
        {
        #pragma omp for
                for (i = 0; i <= end; i++)
                {
                        mchi[i] = pow(10, lmmin + step*i);
                        intOut[i] = intRateDeRoc(r, mchi[i], oper, np, &control_params);
                        // printf("%0.5e\n", intOut[i]);
                        printf("==================================\n");
                        printf("||\t%0.2f %% complete\t||\n", tally * 100. / (end));
                        printf("==================================\n");
                        tally++;
                }
        }

        // char *filename = malloc(strlen("./int_rates/single_op_CE") +  strlen(".dat") + 1);
        // sprintf(filename, "./int_rates/single_op_CE.dat");

        char *filename = "./int_rates/single_op_DR.dat";

        FILE *singleOp = fopen(filename, "w");
        // fprintf(singleOp, "mchi\tOmega\n");
        for(i=0; i<end; i++){
                fprintf(singleOp, "%0.5e\t%0.5e\n", mchi[i], intOut[i]);
        }
        fclose(singleOp);
        // free(filename);

        return 0.;
}

// void method_comps(){
//         int i;
//         int j;
//         int start = -3;
//         int end = 4;
//         int num = 10;
//         int opers = 10;
//         double step = (double)(end - start)/((double) num);

//         double mchi[num+1];
//         double out1[opers+1][num+1];
//         double out2[opers+1][num+1];
//         double out3[opers+1][num+1];

//         // double tally = 0;
//         // double total_evals = 10*((double)num);
//         // #pragma omp parallel
//         // {
//         //         #pragma omp for collapse(2)
//                         for (i=0; i<=num; i++){
//                                 if (opers > 1){
//                                         mchi[i] = pow(10, start + i*step);
//                                         printf("%0.3e\n", mchi[i]);
//                                         printf("-------------------------------------------------------------\n\n");
//                                         for (j = 1; j<=opers; j++){
//                                                 printf("%d\n", j);
//                                                 out1[j-1][i] = crateUappxFull(mchi[i], j, np, &control_params);
//                                                 // printf("M done for %d, %d\n", i, j);
//                                                 // printf("Miser:\t%.5e\n", out1);
//                                                 out2[j-1][i] = cuba_S(mchi[i], j, np, &control_params);
//                                                 // printf("S done for %d, %d\n", i, j);
//                                                 out3[j-1][i] = cuba_D(mchi[i], j, np, &control_params);
//                                                 // printf("D done for %d, %d\n", i, j);
//                                                 printf("=============================================================\n\n");
//                                         }
//                                 }

//                                 else{   
//                                         int oper = 4;
//                                         mchi[i] = pow(10, start + i*step);
//                                         printf("%0.3e\n", mchi[i]);
//                                         printf("-------------------------------------------------------------\n\n");
//                                         out1[0][i] = crateUappxFull(mchi[i], oper, np, &control_params);
//                                         // printf("M done for %d, %d\n", i, 1);
//                                         // printf("Miser:\t%.5e\n", out1);
//                                         out2[0][i] = cuba_S(mchi[i], oper, np, &control_params);
//                                         // printf("S done for %d, %d\n", i, 1);
//                                         out3[0][i] = cuba_D(mchi[i], oper, np, &control_params);
//                                         // printf("D done for %d, %d\n", i, 1);
//                                         printf("=============================================================\n\n");
//                                 }


//                                         // printf("%0.5e\t%0.5e\t%0.5e\t%0.5e\n", mchi[i], out1[j-1][i], out2[j-1][i], out3[j-1][i]);

//                                         // tally++;
//                                         // printf("==================================\n");
//                                         // printf("||\t%0.2f %% complete\t||\n", tally * 100. / (total_evals));
//                                         // printf("==================================\n");
//                                 }
//                 // printf("Ratio M/S:\t%f\n", out1/out2);
//                 // }

//         // printf("\n=============================================================\n\n");

//         // FILE *log =  fopen("./logs/method_comp.dat", "w");
//         for(i = 0; i<=num; i++){
//                 printf("\n%0.5e\t", mchi[i]);
//                 printf("M:\t");
//                 for (j = 0; j<opers; j++){
//                         printf("%0.5e\t", out1[j][i]);
//                 }
//                 printf("\n\t\t\tS:\t");
//                 for (j = 0; j<opers; j++){
//                         printf("%0.5e\t", out2[j][i]);
//                 }
//                 printf("\n\t\t\tD:\t");
//                 for (j = 0; j<opers; j++){
//                         printf("%0.5e\t", out3[j][i]);
//                 }
//                 printf("\n");
//         }
//         //         printf("M:\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n", out1[0][i], out1[1][i], out1[2][i], out1[3][i], out1[4][i], out1[5][i], out1[6][i], out1[7][i], out1[8][i], out1[9][i]);
//         //         printf("S:\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n", out2[0][i], out2[1][i], out2[2][i], out2[3][i], out2[4][i], out2[5][i], out2[6][i], out2[7][i], out2[8][i], out2[9][i]);
//         //         printf("D:\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n", out3[0][i], out3[1][i], out3[2][i], out3[3][i], out3[4][i], out3[5][i], out3[6][i], out3[7][i], out3[8][i], out3[9][i]);

//         //         fprintf(log, "%0.5e\n", mchi[i]);
//         //         fprintf(log, "M:\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n", out1[0][i], out1[1][i], out1[2][i], out1[3][i], out1[4][i], out1[5][i], out1[6][i], out1[7][i], out1[8][i], out1[9][i]);
//         //         fprintf(log, "S:\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n", out2[0][i], out2[1][i], out2[2][i], out2[3][i], out2[4][i], out2[5][i], out2[6][i], out2[7][i], out2[8][i], out2[9][i]);
//         i//         fprintf(log, "D:\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n", out3[0][i], out3[1][i], out3[2][i], out3[3][i], out3[4][i], out3[5][i], out3[6][i], out3[7][i], out3[8][i], out3[9][i]);

//         // }
//         // fclose(log);
//         // printf("=============================================================\n\n");
//         }


//=============================================================================
// SCAN REGIONS
//=============================================================================
void scan_1D_region(double mchi, int oper, int start, int end, float step){

        int num = round((end - start)/step)+1;

        double scan_vals[num];
        double cap_vals[num];

        int i;
        #pragma omp parallel
        {
                #pragma omp for
                        for (i = 0; i<num; i++){

                                double ct = pow(10, start + step * i);

                                struct cont_pars temp_vars = {control_params.z0, control_params.k0, ct};

                                scan_vals[i] = ct;

                                cap_vals[i] = crateFull_T(mchi, kBMeV*1e5, oper, np, &temp_vars);


                                printf("%0.6e\t%0.6e\n", scan_vals[i], cap_vals[i]);
                        }
        }

        FILE *outfile = fopen("../intgration_cut_1D_scan_1_MeV.dat", "w");
        fprintf(outfile, "ct\tC\n");

        for (int i = 0; i<num; i++){
                        for (int j = 0; j<num; j++){
                                fprintf(outfile, "%0.6e\t%0.6e\n", scan_vals[i], cap_vals[i]);
                        }
        }
        fclose(outfile);
}

//=============================================================================

void scan_2D_region(double mchi, int oper, int start, int end, float step){

        int num = round((end - start)/step)+1;

        double z0_vals[num];
        double k0_vals[num];
        double cap_vals[num][num];

        #pragma omp parallel
        {
                #pragma omp for collapse(2)
                        for (int i = 0; i<num; i++){
                                for (int j = 0; j<num; j++){
                                        double k0 = pow(10, start + step * i);
                                        double z0 = pow(10, start + step * j);

                                        struct cont_pars temp_vars = {z0, k0, 100};

                                        cap_vals[i][j] = crateFull(mchi, oper, np, &temp_vars);

                                        k0_vals[i] = k0;
                                        z0_vals[j] = z0;

                                        printf("%0.6e\t%0.6e\t%0.6e\n", k0_vals[i], z0_vals[j], cap_vals[i][j]);
                                }
                        }
        }

        FILE *outfile = fopen("../intgration_cut_scan_1_MeV.dat", "w");
        fprintf(outfile, "k0\tz0\tC\n");

        for (int i = 0; i<num; i++){
                        for (int j = 0; j<num; j++){
                                fprintf(outfile, "%0.6e\t%0.6e\t%0.6e\n", k0_vals[i], z0_vals[j], cap_vals[i][j]);
                        }
        }
        fclose(outfile);
}

//=============================================================================
// CAPTURE ROUTINES
//=============================================================================

double crate_complete(double lmmin, double lmmax, double step)
{
        int i1, i2;                                       // mass loop
        int j;                                            // operator loop

        float lmmin1 = lmmin;
        float lmmax1 = log10(6e4);

        float lmmin2 = lmmax1 + step;
        float lmmax2 = lmmax;

        int end1 = round((lmmax1 - lmmin1) / step) + 1; //number of mass points
        int end2 = round((lmmax2 - lmmin2) / step) + 1; //number of mass points
        int end_tot = end1+end2;
        double capOut[11][end_tot];

        double mchi[end_tot];

        printf("Computing capture rate over full mass range from 1e%0.f MeV to 1e%0.f MeV\n", lmmin, lmmax);

        int tally = 0;
        double total_evals = (end_tot) * 10 - 1;

        #pragma omp parallel
                {
                #pragma omp for collapse(2)
                        for (i1 = 0; i1 < end1; i1++)
                        {
                                for (j = 1; j <= 10; j++)
                                {
                                        mchi[i1] = pow(10, lmmin1 + step * i1);
                                        capOut[j - 1][i1] = crateUappxFull(mchi[i1], j, np, &control_params);
                                        printf("==================================\n");
                                        printf("||\t%0.2f %% complete\t||\n", tally * 100. / (total_evals));
                                        printf("==================================\n");
                                        tally++;
                                }
                        }
                j = 0;
                #pragma omp for collapse(2)
                        for (i2 = 0; i2 < end2; i2++)
                        {
                                for (j = 1; j <= 10; j++)
                                {
                                        mchi[i2 + end1] = pow(10, lmmin2 + step * i2);
                                        capOut[j - 1][i2 + end1] = crateMultiScatt(mchi[i2 + end1], j, np);
                                        // capOut[j - 1][i2 + end1] = crateHighMass(mchi[i2 + end1], j, np);
                                        // printf("m = %0.5e,  exp = %0.5e,  i = %d\n", mchi[i2+end1], lmmin2 + step * i2, i2);
                                        printf("==================================\n");
                                        printf("||\t%0.2f %% complete\t||\n", tally * 100. / (total_evals));
                                        printf("==================================\n");
                                        tally++;
                                }
                        }
                }

        char *filename = malloc(strlen("./crate_all/complete///crate_") + strlen(label) + strlen(set_WD) + strlen(elem) + strlen(".dat") + 1);
        sprintf(filename, "./crate_all/complete/%s/%s/crate_%s.dat", label, elem, set_WD);

        FILE *capoutfile = fopen(filename, "w");
        fprintf(capoutfile, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                "mDM_(MeV)", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10");
        int i;
        for (i = 0; i < end_tot; i++)
        {
                fprintf(capoutfile, "%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n",
                        mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);

                // printf("%0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e\n",
                //        mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);
                // printf("m = %0.5e,  i = %d\n", mchi[i], i);
        }

        fclose(capoutfile);
        free(filename);

        return 0.;
}

//=============================================================================
// CAPTURE RATES
//=============================================================================


double crate_singleOp_Full(int oper, double lmmin, double lmmax, double step)
{

        int i; // mass loop
        int end = round((lmmax - lmmin) / step)+1; //number of mass points
        double capOut[end];
        double mchi[end];
	int tally = 0;
        #pragma omp parallel
        {
        #pragma omp for
                for (i = 0; i <= end; i++)
                {
                        mchi[i] = pow(10, lmmin + step*i);
                        capOut[i] = crateUappxFull(mchi[i], oper, np, &control_params);
                        printf("==================================\n");
                        printf("||\t%0.2f %% complete\t||\n", tally * 100. / (end));
                        printf("==================================\n");
                        tally++;
                }
        }

        char *filename = malloc(strlen("./crate_single/single_op_full_") + strlen(set_WD) + strlen(elem) +  strlen("__.dat") + 1);
        sprintf(filename, "./crate_single/single_op_full_%s_%s.dat", set_WD, elem);

        FILE *singleOp = fopen(filename, "w");
        for(i=0; i<end; i++){
                fprintf(singleOp, "%0.5e\t%0.5e\n", mchi[i], capOut[i]);
        }
        fclose(singleOp);
        free(filename);

        return 0.;
}

double crate_singleOp_NoPB(int oper, double lmmin, double lmmax, double step)
{


        int i; // mass loop
        int end = round((lmmax - lmmin) / step) + 1; //number of mass points
        double capOut[end];
        double mchi[end];

        int tally = 0;

        #pragma omp parallel
                {
        #pragma omp for
                        for (i = 0; i <= end; i++)
                        {
                                mchi[i] = pow(10, lmmin + step * i);
                                capOut[i] = crateNoPb(mchi[i], oper, np);
                                printf("==================================\n");
                                printf("||\t%0.2f %% complete\t||\n", tally * 100. / (end));
                                printf("==================================\n");
                                tally++;
                        }
                }

        char *filename = malloc(strlen("./crate_single/single_op_crateNoPB_") + strlen(set_WD) + strlen(elem) +  strlen("__.dat") + 1);
        sprintf(filename, "./crate_single/single_op_crateNoPB_%s_%s_.dat", set_WD, elem);

        FILE *singleOp = fopen(filename, "w");
        for (i = 0; i < end; i++)
        {
                fprintf(singleOp, "%0.5e\t%0.5e\n", mchi[i], capOut[i]);
        }
        fclose(singleOp);
        free(filename);

        return 0.;
}

double crate_singleOp_CollEff(int oper, double lmmin, double lmmax, double step)
{

        int i; // mass loop
        int end = round((lmmax - lmmin) / step)+1; //number of mass points
        double capOut[end];
        double mchi[end];
	int tally = 0;
        #pragma omp parallel
        {
        #pragma omp for
                for (i = 0; i <= end; i++)
                {
                        mchi[i] = pow(10, lmmin + step*i);
                        capOut[i] = crateCollEff(mchi[i], oper, np, &control_params);
                        printf("==================================\n");
                        printf("||\t%0.2f %% complete\t||\n", tally * 100. / (end));
                        printf("==================================\n");
                        tally++;
                }
        }

        char *filename = malloc(strlen("./crate_single/single_op_coll_eff_") + strlen(set_WD) + strlen(elem) +  strlen("__.dat") + 1);
        sprintf(filename, "./crate_single/single_op_coll_eff_%s_%s.dat", set_WD, elem);

        FILE *singleOp = fopen(filename, "w");
        for(i=0; i<end; i++){
                fprintf(singleOp, "%0.5e\t%0.5e\n", mchi[i], capOut[i]);
        }
        fclose(singleOp);
        free(filename);

        return 0.;
}


double crate_singleOp_finiteT_full(int oper, double lmmin, double lmmax, double temp, double step)
{


        int i; // mass loop
        int end = round((lmmax - lmmin) / step) + 1; //number of mass points
        double capOut[end];
        double mchi[end];
        int tally = 0;

        #pragma omp parallel
        {
        #pragma omp for
                for (i = 0; i <= end; i++)
                {
                        mchi[i] = pow(10, lmmin + step * i);
                        capOut[i] = crateFull_T(mchi[i], temp, oper, np, &control_params);
                        printf("==================================\n");
                        printf("||\t%0.2f %% complete\t||\n", tally * 100. / (end));
                        printf("==================================\n");
                        tally++;
                }
        }

        char *filename = malloc(strlen("./crate_single/single_op_finite_T_full_") + strlen(set_WD) + strlen(elem) +  strlen("__.dat") + 1);
        sprintf(filename, "./crate_single/single_op_finite_T_full_%s_%s_.dat", set_WD, elem);

        FILE *singleOp = fopen(filename, "w");
        for (i = 0; i < end; i++)
        {
                fprintf(singleOp, "%0.5e\t%0.5e\n", mchi[i], capOut[i]);
        }
        fclose(singleOp);
        free(filename);

        return 0.;
}

double crate_all_Full(double lmmin, double lmmax, double step)
{

        int i; // mass loop
        int end = round((lmmax - lmmin) / step)+1;  //number of mass points
        int j; // operator loop
        double capOut[11][end];

        double mchi[end];

        int tally = 0;
        double total_evals = end*11;



        #pragma omp parallel
        {
                #pragma omp for collapse(2)
                for (i = 0; i < end; i++)
                {
                        for (j = 1; j <= 10; j++)
                        {
                                mchi[i] = pow(10, lmmin + step*i);
                                capOut[j - 1][i] = crateFull(mchi[i], j, np, &control_params);
                                printf("==================================\n");
                                printf("||\t%0.2f %% complete\t||\n", tally * 100. / (total_evals));
                                printf("==================================\n");
                                tally++;
                        }
                }
        }

        // char *filename = malloc(strlen("./crate_all/crate_operdn_full_") + strlen(set_WD) + strlen(elem) +  strlen(".dat") + 1);
        char *filename = malloc(strlen("./crate_all/full/pauli///crate_") + strlen(label) + strlen(set_WD) + strlen(elem) +  strlen(".dat") + 1);
        sprintf(filename, "./crate_all/full/pauli/%s/%s/crate_%s.dat", label, elem, set_WD);

        FILE *capoutfile = fopen(filename, "w");
        fprintf(capoutfile, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                "mDM_(MeV)", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10");

        for (i = 0; i < end; i++)
        {
                fprintf(capoutfile, "%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n",
                        mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);

                printf("%0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e\n",
                mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);
        }

        fclose(capoutfile);
        free(filename);
        return 0.;
}

double crate_all_NoPB(double lmmin, double lmmax, double step)
{
        int i; // mass loop
        int end = round((lmmax - lmmin) / step)+1; //number of mass points
        int j; // operator loop
        double capOut[11][end];

        double mchi[end];
        int tally = 0;
        double total_evals = end*11;



        #pragma omp parallel
        {
                #pragma omp for collapse(2)
                for (i = 0; i < end; i++)
                {
                        for (j = 1; j <= 10; j++)
                        {
                                mchi[i] = pow(10, lmmin + step*i);
                                capOut[j - 1][i] = crateNoPb(mchi[i], j, np);
                                printf("==================================\n");
                                printf("||\t%0.2f %% complete\t||\n", tally * 100. / (total_evals));
                                printf("==================================\n");
                                tally++;
                        }
                }
        }

        char *filename = malloc(strlen("./crate_all/full/no_pauli///crate_") + strlen(label) + strlen(set_WD) + strlen(elem) +  strlen(".dat") + 1);
        sprintf(filename, "./crate_all/full/no_pauli/%s/%s/crate_%s.dat", label, elem, set_WD);

        FILE *capoutfile = fopen(filename, "w");
        fprintf(capoutfile, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                "mDM_(MeV)", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10");

        for (i = 0; i < end; i++)
        {
                fprintf(capoutfile, "%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n",
                                mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);

                printf("%0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e\n",
                        mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);
        }

        fclose(capoutfile);
        free(filename);

        return 0.;

}

//=============================================================================
// CAPTURE RATES U APPX
//=============================================================================

double crate_singleOp_Uappx_Full(int oper, double lmmin, double lmmax, double step)
{


        int i; // mass loop
        int end = round((lmmax - lmmin) / step) + 1; //number of mass points
        double capOut[end];
        double mchi[end];
        int tally = 0;
        #pragma omp parallel
        {

        #pragma omp for
                for (i = 0; i <= end; i++)
                {
                        mchi[i] = pow(10, lmmin + step * i);
                        capOut[i] = crateUappxFull(mchi[i], oper, np, &control_params);
                        printf("==================================\n");
                        printf("||\t%0.2f %% complete\t||\n", tally * 100. / (end));
                        printf("==================================\n");
                        tally++;
                }
        }

        char *filename = malloc(strlen("./crate_single/single_op_crateuappx_") + strlen(set_WD) + strlen(elem) +  strlen("__.dat") + 1);
        sprintf(filename, "./crate_single/single_op_crateuappx_%s_%s.dat", set_WD, elem);
 
        FILE *singleOp = fopen(filename, "w");
        for (i = 0; i < end; i++)
        {
                fprintf(singleOp, "%0.5e\t%0.5e\n", mchi[i], capOut[i]);
        }
        fclose(singleOp);
        free(filename);

        return 0.;
}

double crate_all_Uappx_Full(double lmmin, double lmmax, double step)
{

        int i; // mass loop
        int end = round((lmmax - lmmin) / step) + 1; //number of mass points
        int j; // operator loop
        int tally=0;
        int total_evals = 11 * end;
        double capOut[11][end];

        double mchi[end];



        #pragma omp parallel
        {
        #pragma omp for collapse(2)
                for (i = 0; i < end; i++)
                {
                        for (j = 1; j <= 10; j++)
                        {
                                mchi[i] = pow(10, lmmin + step * i);
                                capOut[j - 1][i] = crateUappxFull(mchi[i], j, np, &control_params);
                                printf("m_chi = %0.3e, oper = D%d, C = %0.5e\n", mchi[i], j, capOut[j-1][i]);
                                printf("==================================\n");
                                printf("||\t%0.2f %% complete\t||\n", tally * 100. / (total_evals));
                                printf("==================================\n");
                                tally++;
                        }
                }
        }

        char *filename = malloc(strlen("./crate_all/Uappx/pauli///crate_") + strlen(label) + strlen(set_WD) + strlen(elem) +  strlen(".dat") + 1);
        sprintf(filename, "./crate_all/Uappx/pauli/%s/%s/crate_%s.dat", label, elem, set_WD);

        FILE *capoutfile = fopen(filename, "w");
        fprintf(capoutfile, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                "mDM_(MeV)", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10");

        for (i = 0; i < end; i++)
        {
                fprintf(capoutfile, "%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n",
                        mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);

                printf("%0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e\n",
                       mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);
        }

        fclose(capoutfile);
        free(filename);

        return 0.;
}

double crate_all_Uappx_NoPB(double lmmin, double lmmax, double step)
{
        int i; // mass loop
        int end = round((lmmax - lmmin) / step) + 1; //number of mass points
        int j; // operator loop
        double capOut[11][end];
        int tally = 0;
        int total_evals = 11*end;

        double mchi[end];



        #pragma omp parallel
        {
        #pragma omp for collapse(2)
                for (i = 0; i < end; i++)
                {
                        for (j = 1; j <= 10; j++)
                        {
                                mchi[i] = pow(10, lmmin + step * i);
                                capOut[j - 1][i] = crateUappxNoPb(mchi[i], j, np);

                                printf("==================================\n");
                                printf("||\t%0.2f %% complete\t||\n", tally * 100. / (total_evals));
                                printf("==================================\n");
                                tally++;
                        }
                }
        }

        char *filename = malloc(strlen("./crate_all/Uappx/no_pauli///crate_") + strlen(label) + strlen(set_WD) + strlen(elem) +  strlen("__.dat") + 1);
        sprintf(filename, "./crate_all/Uappx/no_pauli/%s/%s/crate_%s.dat", label, elem, set_WD);

        FILE *capoutfile = fopen(filename, "w");
        fprintf(capoutfile, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                "mDM_(MeV)", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10");

        for (i = 0; i < end; i++)
        {
                fprintf(capoutfile, "%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n",
                        mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);

                printf("%0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e\n",
                       mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);
        }

        fclose(capoutfile);
        free(filename);
        return 0.;

        return 0.;
}

//=====================================================
// CAPTURE RADIAL PROFILES
//=====================================================

// double cap_prof_full(double rmin, double rmax, double mchi, int num_pts)
// {
//         int opers[2] = {3, 5};
//         double radcapOut[2][7][num_pts];
//         double rad[num_pts];
//         int i, j;
//         int tally = 0;

//         int num_mass = 7;

//         double mass_vals[7] = {10, 50, 100, 500, 1000, 5000, 10000};
//         char *mass_names[7] = {"10", "50", "100", "500", "1000", "5000", "10000"};

// #pragma omp parallel
//         {
//         #pragma omp for collapse(3)
//                 for (int m = 0; m < num_mass; m++)
//                 {
//                         for (i = 0; i < num_pts; i++)
//                         {
//                                 for ( j = 0; j <= 1; j++ )
//                                 {
//                                         rad[i] = rmin + i*(rmax - rmin)/(num_pts - 1);
//                                         radcapOut[j][m][i] = crateRaProfFull(rad[i], mass_vals[m], opers[j], np);

//                                         printf("==================================\n");
//                                         printf("||\t%0.2f %% complete\t||\n", tally * 100. /10./ (num_pts));
//                                         printf("==================================\n");
//                                         tally++;
//                                 }
//                         }
//                 }
//         }

//         for (int m = 0; m < num_mass; m++){
//         char *filename = malloc(strlen("./radial_profiles/crad_prof_full_") + strlen(set_WD) + strlen("_") + strlen(mass_names[m]) + strlen(".dat") + 1);
//         sprintf(filename, "./radial_profiles/crad_prof_full_%s_%s.dat", set_WD, mass_names[m]);

//         FILE *radprofout = fopen(filename, "w");

//         fprintf(radprofout, "%s\t%s\n",
//                 "r/R", "D3", "D5");

//         for (i = 0; i < num_pts; i++){
//                 fprintf(radprofout, "%0.5e\t%0.5e\t%0.5e\n",
//                  rad[i], radcapOut[0][m][i], radcapOut[1][m][i]);
//         }

//         // fprintf(radprofout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
//         //         "r/R", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10");

//         // for (i = 0; i <= num_pts; i++){
//         //         fprintf(radprofout, "%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n",
//         //          rad[i], radcapOut[0][i], radcapOut[1][i], radcapOut[2][i], radcapOut[3][i], radcapOut[4][i], radcapOut[5][i], radcapOut[6][i], radcapOut[7][i], radcapOut[8][i], radcapOut[9][i]);
//         // }

//         fclose(radprofout);
//         free(filename);
//         }
//         return 0.;
// }

// double cap_prof_nopb(int oper, double rmin, double rmax, double mchi, double step)
// {


//         int end = round((rmax - rmin) / step) + 1; //number of mass points
//         double radcapOut[end];
//         double rad[end];
//         int i;
//         printf("%d\n", end);

// #pragma omp parallel
//         {
// #pragma omp for
//         for (i = 0; i < end; i++)
//         {
//                 rad[i] = rmin + i*step;
//                 radcapOut[i] = crateRaProfNoPb(rad[i], mchi, oper, np);
//         }
//         }

//         FILE *radprofout = fopen("./radial_profiles/crad_prof_nopb.dat", "w");

//         for (i=0; i< end; i++){
//                 fprintf(radprofout, "%0.5e\t%0.5e\n", rad[i], radcapOut[i]);
//                 printf("%0.5e\t%0.5e\n", rad[i], radcapOut[i]);
//         }

//         fclose(radprofout);

//         return 0.;
// }

// double cap_prof_appx(int oper, double rmin, double rmax, double mchi, double step)
// {


//         int end = round((rmax - rmin) / step) + 1; //number of mass points
//         double radcapOut[end];
//         double rad[end];
//         int i;
//         int iters = 0;

// #pragma omp parallel
//         {
// #pragma omp for
//         for (i = 0; i < end; i++)
//         {
//                 rad[i] = rmin + i*step;
//                 radcapOut[i] = crateRaProfAppx(rad[i], mchi, oper, np, iters);
//         }
//         }

//         FILE *radprofout = fopen("./radial_profiles/crad_prof_appx.dat", "w");

//         for (i=0; i< end; i++){
//                 fprintf(radprofout, "%0.5e\t%0.5e\n", rad[i], radcapOut[i]);
//                 printf("%0.5e\t%0.5e\n", rad[i], radcapOut[i]);
//         }

//         fclose(radprofout);

//         return 0.;
// }

//=====================================================
// FINITE TEMPERATURE CAPTURE
//=====================================================

double crate_all_finiteT_NoPB(double lmmin, double lmmax, double temp, double step)
{
        int i; // mass loop
        int end = round((lmmax - lmmin) / step) + 1; //number of mass points
        int j; // operator loop
        double capOut[11][end];

        double mchi[end];
        int tally = 0;
        double total_evals = end * 11;



#pragma omp parallel
        {
#pragma omp for collapse(2)
                for (i = 0; i < end; i++)
                {
                        for (j = 1; j <= 10; j++)
                        {
                                mchi[i] = pow(10, lmmin + step * i);
                                capOut[j - 1][i] = crateNoPb_T(mchi[i], temp, j, np, &control_params);


                                printf("==================================\n");
                                printf("||\t%0.2f %% complete\t||\n", tally * 100. / (total_evals));
                                printf("==================================\n");
                                tally++;
                        }
                }
        }

        char *filename = malloc(strlen("./crate_finite_temp/full/no_pauli///crate_") + strlen(label) + strlen(set_WD) + strlen(elem) +  strlen(".dat") + 1);
        sprintf(filename, "./crate_finite_temp/full/no_pauli/%s/%s/crate_%s.dat",label, elem, set_WD);

        FILE *capoutfile = fopen(filename, "w");
        fprintf(capoutfile, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                "mDM_(MeV)", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10");

        for (i = 0; i < end; i++)
        {
                fprintf(capoutfile, "%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n",
                        mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);

                printf("%0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e\n",
                       mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);
        }

        fclose(capoutfile);
        free(filename);

        return 0.;
}

double crate_all_finiteT_Full(double lmmin, double lmmax, double temp, double step)
{

        int i; // mass loop
        int end = round((lmmax - lmmin) / step) + 1; //number of mass points
        int j; // operator loop
        double capOut[11][end];

        double mchi[end];

        int tally = 0;
        double total_evals = end * 11;

        printf("Computing Finite Temperature Full Calculation...\n");


#pragma omp parallel
        {
#pragma omp for collapse(2)
                for (i = 0; i < end; i++)
                {
                        for (j = 1; j <= 10; j++)
                        {
                                mchi[i] = pow(10, lmmin + step * i);
                                capOut[j - 1][i] = crateFull_T(mchi[i], temp, j, np, &control_params);


                                printf("==================================\n");
                                printf("||\t%0.2f %% complete\t||\n", tally * 100. / (total_evals));
                                printf("==================================\n");
                                tally++;
                        }
                }
        }

        char *filename = malloc(strlen("./crate_finite_temp/full///pauli/crate_") + strlen(label) + strlen(set_WD) + strlen(elem) +  strlen(".dat") + 1);
        sprintf(filename, "./crate_finite_temp/full/pauli/%s/%s/crate_%s.dat", label, elem, set_WD);

        FILE *capoutfile = fopen(filename, "w");
        fprintf(capoutfile, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                "mDM_(MeV)", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10");

        for (i = 0; i < end; i++)
        {
                fprintf(capoutfile, "%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n",
                        mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);

                printf("%0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e\n",
                       mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);
        }

        fclose(capoutfile);
        free(filename);

        return 0.;
}

double crate_all_finiteT_Uappx(double lmmin, double lmmax, double temp, double step)
{

        int i;                                       // mass loop
        int end = round((lmmax - lmmin) / step) + 1; //number of mass points
        int j;                                       // operator loop
        double capOut[11][end];

        double mchi[end];

        int tally = 0;
        double total_evals = end * 11;

        printf("Computing Finite Temperature PB with uchi approximation...\n");

#pragma omp parallel
        {
#pragma omp for collapse(2)
                for (i = 0; i < end; i++)
                {
                        for (j = 1; j <= 10; j++)
                        {
                                mchi[i] = pow(10, lmmin + step * i);
                                capOut[j - 1][i] = crateUappx_T(mchi[i], temp, j, np, &control_params);


                                printf("==================================\n");
                                printf("||\t%0.2f %% complete\t||\n", tally * 100. / (total_evals));
                                printf("==================================\n");
                                tally++;
                        }
                }
        }

        char *filename = malloc(strlen("./crate_finite_temp/Uappx/pauli///crate_") + strlen(label) + strlen(set_WD) + strlen(elem) +  strlen(".dat") + 1);
        sprintf(filename, "./crate_finite_temp/Uappx/pauli/%s/%s/crate_%s.dat", label, elem, set_WD);

        FILE *capoutfile = fopen(filename, "w");
        fprintf(capoutfile, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                "mDM_(MeV)", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10");

        for (i = 0; i < end; i++)
        {
                fprintf(capoutfile, "%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n",
                        mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);

                printf("%0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e  %0.5e\n",
                       mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);
        }

        fclose(capoutfile);
        free(filename);

        return 0.;
}
//=====================================================
// MULTISCATTER REGIEME
//=====================================================

double crate_all_multiscat(double lmmin, double lmmax, double step)
{

        int i; // mass loop
        int end = round((lmmax - lmmin) / step)+1; //number of mass points
        int j; // operator loop
        double capOut[11][end];

        double mchi[end];

        int tally = 0;
        double total_evals = end*11;



        #pragma omp parallel
        {
                #pragma omp for collapse(2)
                for (i = 0; i < end; i++)
                {
                        for (j = 1; j <= 10; j++)
                        {
                                mchi[i] = pow(10, lmmin + step*i);
                                capOut[j - 1][i] = crateMultiScatt(mchi[i], j, np);
                                // capOut[j - 1][i] = crateHighMass(mchi[i], j, np);


                                printf("==================================\n");
                                printf("||\t%0.2f %% complete\t||\n", tally * 100. / (total_evals));
                                printf("==================================\n");
                                tally++;
                        }
                }
        }

        char *filename = malloc(strlen("./crate_all/multiscatt///crate_") + strlen(label) + strlen(set_WD) + strlen(elem) +  strlen(".dat") + 1);
        sprintf(filename, "./crate_all/multiscatt/%s/%s/crate_%s.dat", label, elem, set_WD);

        FILE *capoutfile = fopen(filename, "w");
        fprintf(capoutfile, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                "mDM_(MeV)", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10");

        for (i = 0; i < end; i++)
        {
                fprintf(capoutfile, "%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\n",
                        mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);

                printf("%0.3e  %0.3e  %0.3e  %0.3e  %0.3e  %0.3e  %0.3e  %0.3e  %0.3e  %0.3e  %0.3e\n",
                mchi[i], capOut[0][i], capOut[1][i], capOut[2][i], capOut[3][i], capOut[4][i], capOut[5][i], capOut[6][i], capOut[7][i], capOut[8][i], capOut[9][i]);
        }

        fclose(capoutfile);
        free(filename);

        return 0.;
}
