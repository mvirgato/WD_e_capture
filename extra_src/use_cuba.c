#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "globals.h"
#include "use_cuba.h"
#include "cuba.h"
#include "crate_apprx_funcs.h"
#include "cap_interp.h"

//===================================================================================

#define NDIM 4
#define NCOMP 1
#define NVEC 1
#define EPSREL 1e-2
#define EPSABS 1e-12
#define VERBOSE 0
#define LAST 4
#define SEED 0
#define MINEVAL 10000
#define MAXEVAL 500000

#define NSTART 10000
#define NINCREASE 1000
#define NBATCH 10000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 10000
#define NMIN 500
#define FLATNESS 50.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

//===================================================================================

int cuba_int_uappx(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {

  struct crate_params *params = (struct crate_params *)userdata;

  double mdm = (params->mchi);
  int np = (params->npts);
  int oper = (params->oper);
  double k0 = (params->k0);
  double z0 = (params->z0);

  #define r xx[0]
  #define E xx[1]
  #define s xx[2]
  #define t xx[3]

  #define f ff[0]

  double muf = muFeinterp(r, np);
  double ne  = nEinterp(r, np);
  double B   = Binterp(r, np);

  f = 1e30*integrandUappx(r, E, s, t, muf, ne, B, mdm, oper, k0, z0);
  // f = r*E*s*t;
  // printf("%0.5e\n",f);

  return 0;
}

//===================================================================================

double cuba_S(double mchi, int oper, int npts, void *cont_vars) {

    struct cont_pars *contp  = (struct cont_pars *)cont_vars;

    double z0 = (contp->z0);
    double k0 = (contp->k0);

    struct crate_params USERDATA = {mchi, oper, npts, z0, k0};

    int comp, nregions, fail;
    long long neval;
    cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
    
    // printf("\n-------------------- Suave test --------------------\n");
    
    llSuave(NDIM, NCOMP, cuba_int_uappx, &USERDATA, NVEC,
        EPSREL, EPSABS, VERBOSE | LAST, SEED,
        MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
        STATEFILE, SPIN,
        &nregions, &neval, &fail, integral, error, prob);
    
    // printf("SUAVE RESULT:\tnregions %d\tneval %lld\tfail %d\n",
    //     nregions, neval, fail);
    for( comp = 0; comp < NCOMP; ++comp ){
        printf("SUAVE RESULT:\t%.5e +- %.3e\tp = %.3e\n", (double)integral[comp]*CoeffUappx(mchi), (double)error[comp]*CoeffUappx(mchi), (double)prob[comp]);
    }
    // printf("Full cap rate is:\t%0.5e\n", (double)integral[0]*CoeffUappx(mchi));
    // printf("Suave:\t%0.5e\n", (double)integral[0]*CoeffUappx(mchi));
    


  return (double)integral[0]*CoeffUappx(mchi);
}

double cuba_D(double mchi, int oper, int npts, void *cont_vars)
{

  struct cont_pars *contp = (struct cont_pars *)cont_vars;

  double z0 = (contp->z0);
  double k0 = (contp->k0);

  struct crate_params USERDATA = {mchi, oper, npts, z0, k0};

  int comp, nregions, fail;
  long long int neval;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  // printf("\n------------------- Divonne test -------------------\n");

  llDivonne(NDIM, NCOMP, cuba_int_uappx, &USERDATA, NVEC,
          EPSREL, EPSABS, VERBOSE, SEED,
          MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
          BORDER, MAXCHISQ, MINDEVIATION,
          NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
          STATEFILE, SPIN,
          &nregions, &neval, &fail, integral, error, prob);

  // printf("DIVONNE RESULT:\tnregions %d\tneval %lld\tfail %d\n",
  //        nregions, neval, fail);
  for (comp = 0; comp < NCOMP; ++comp){
    printf("DIVONNE RESULT:\t%.5e +- %.3e\tp = %.3e\n", (double)integral[comp]*CoeffUappx(mchi), (double)error[comp]*CoeffUappx(mchi), (double)prob[comp]);
           }
  // printf("Suave:\t%0.5e\n", (double)integral[0]*CoeffUappx(mchi));

  return (double)integral[0] * CoeffUappx(mchi);
}

double cuba_C(double mchi, int oper, int npts, void *cont_vars)
{

  struct cont_pars *contp = (struct cont_pars *)cont_vars;

  double z0 = (contp->z0);
  double k0 = (contp->k0);

  struct crate_params USERDATA = {mchi, oper, npts, z0, k0};

  int comp, nregions, fail;
  long long int neval;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

#if 1
  printf("\n-------------------- Cuhre test --------------------\n");

  llCuhre(NDIM, NCOMP, cuba_int_uappx, &USERDATA, NVEC,
        EPSREL, EPSABS, VERBOSE | LAST,
        MINEVAL, MAXEVAL, KEY,
        STATEFILE, SPIN,
        &nregions, &neval, &fail, integral, error, prob);

  // printf("CUHRE RESULT:\tnregions %d\tneval %lld\tfail %d\n",
  //        nregions, neval, fail);
  for (comp = 0; comp < NCOMP; ++comp)
    printf("CUHRE RESULT:\t%.5e +/- %.5e\tp = %.3e\n",
           (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif
  // printf("Suave:\t%0.5e\n", (double)integral[0]*CoeffUappx(mchi));

  return (double)integral[0] * CoeffUappx(mchi);
}