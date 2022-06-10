#ifndef GLOBALS_H
#define GLOBALS_H

static double sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)

static double cubarg;
#define CUB(a) ((cubarg = (a)) == 0.0 ? 0.0 : cubarg * cubarg * cubarg)

double mt;

//unit conversions
static const double mTOcm = 100.;                           // m to cm
static const double cmTOm = 1. / 100.;                      // cm to m
static const double kmTOm = 1000.;                          // km to m
static const double mTOkm = 1. / 1000.;                     // m to km
static const double mTOpm = 1e12;                           // m to pm
static const double pmTOm = 1e-12;                          // pm to m
static const double mTOfm = 1e15;                           // m to fm
static const double fmTOm = 1e-15;                          // fm to m
static const double kgTOg = 1000.;                          // kg to g
static const double gTOkg = 1. / 1000.;                     // g to kg
static const double PaTOdynpcm2 = 10.;                      // (*Pascals to dyne/cm^2*)
static const double inveVTOm = 197.3269788 * fmTOm;         // (*MeV m*)
static const double mTOinveV = 1. / (197.3269788 * fmTOm);  // (*(MeV m)^-1*)
static const double inveVTOs = 6.582119514e-22;
static const double sTOinveV = 1. / (6.582119514e-22);

// general constants 
static const double cspeed  = 299792458.;       // speed of light m/s
static const double hplanck = 6.62607015e-34;   // plancks cosnst. J s
static const double vev     = 246000.;          // Higgs vev MeV
static const double Msol    = 1.989e30;         // Mass of Sun kg
static const double GNewt   = 6.67430e-11;      // Newtons constant SI
static const double pi      = M_PI;             // Pi
static const double kBMeV   = 8.61733034e-11;   // Boltzmann's const. in MeV/K
static const double hbar    = 6.582119569e-22;  // reduced plancks const. MeV s
static const double echarge = 0.302862;         // electric charge in natural units

// halo parameters
static const double vs      = 20.; // star velocity km/s
static const double vd      = 8.; // DM dispersion velocity km/s
static const double rhoDM   = 798e3; // DM density in neighborhood  MeV / cm^3

// target masses
static const double mN = 939; // neutron mass (*MeV*)
static const double mP = 938; // proton mass (*MeV*)
static const double mE = 0.511; // electron mass (*MeV*)
static const double mmu = 105.66; // muon mass (*MeV*)

// Capture rate integration structures
struct cont_pars {double z0; double k0; double ct;};
struct crate_params {double mchi; int oper; int npts; double z0; double k0;};
struct crate_params_T {double mchi; double temp; int oper; int npts; double z0; double k0; double ct;};
struct radprof_params {double rad; double mchi; int oper; int npts;};

// Control Variables
// static const double k0 = 0.0005;
// static const double z0 = 0.05;

#endif