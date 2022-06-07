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


double q0_tr(double uchi, double mchi, double B, double muFe, double soln){
    // Energy transfer

    double E1 = mchi*(1 + uchi*uchi/2.)/sqrt(B); // Initial DM energy

    double _q0 = E1 - soln - 0.5 * mchi * uchi * uchi / sqrt(B);

    return _q0;
}

double q_tr_2(double q0, double t){

    // Squared 3-momentum transfer

    double _q2 = q0*q0 - t; // squared momentum transfer

    return _q2;

}

double Eplus(double t, double q0, double q, double muFe){
    // Kinematic factor (B.8) in Coll. Effects paper

    return fmax(muFe, 0.5 * ( q0 + q * sqrt(1.0 - 4.0 * mt*mt / t)));
}

