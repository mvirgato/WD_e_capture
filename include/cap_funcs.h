void display_results (char *title, double result, double error);
void system_setup(double mT, char *WDname, char *eos_type, char *elem, char *vstar, int *np);
double *logspace(double a, double b, int n, double u[]);

// cap funcs
double fMB(double uchi);
double beta_cap(double s, double mchi);
double gamma_cap(double s, double mchi); // gamma in cap funcs
double solnFunc(double Eu, double s, double t, double uchi, double mchi, double B, double muFe);
double HeavisidePB(double Eu, double s, double t, double uchi, double mchi, double B, double muFe, double soln);
double HeavisideUchi(double Eu, double s, double t, double uchi, double mchi, double B, double muFe, double soln);
double HeavisidePhaseSpace(double Eu, double s, double B, double _smin, double _smax);
double nfree(double muFe);
double zeta(double nE, double muFe);
double Coeff(double mchi);


// Integration limits
double collcor(double mchi, double z0);
double tmin(double s, double mchi);
double tmax();
double smin(double Eu, double uchi, double B, double mchi);
double smax(double Eu, double uchi, double B, double mchi);
double Eumin(double muFe, double mchi, double k0);
double Eumax(double muFe);
double umax();
double rmin(double mchi);
double rmax();

// Integrands

double integrandNoPB(double r, double uchi, double Eu, double sp, double tp, double muFe, double nE, double B, double mchi, int oper);
double integrandFull(double r, double uchi, double Ep, double sp, double tp, double muFe, double nE, double B, double mchi, int oper, double z0, double k0);
// double monteIntegrandNoPB(double x[], size_t dim, void * p);
// double monteIntegrandFull(double x[], size_t dim, void *p);

// Capture radial profile integrands 
// double monteIntegrandCrateRadFull(double x[], size_t dim, void *p);
// double monteIntegrandCrateRadNoPB(double x[], size_t dim, void *p);

// Integrators
double crateNoPb(double mchi, int oper, int npts);
double crateFull(double mchi, int oper, int npts, void *cont_vars);

// Cross sections
double dXS_const(double s, double t, double mchi);
double opersdXS(double s, double t, double mchi, int num);


