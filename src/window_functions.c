#include <math.h>
#include <gsl/gsl_spline.h>
#include "macros_and_structures.h"
#include "functions.h"



//////////////////////////////////
// GLOBAL VARIABLES FROM MAIN.C //
//////////////////////////////////
/* Struct Global Variables */
extern struct globalVariables GV;
extern gsl_interp_accel *acc;
extern gsl_spline *splineRe;
extern gsl_spline *splineIm;
extern gsl_spline *splineMag2;

double sinc(double x)
{
  if(fabs(x)<1.0e-10)
    return 1.0;
  else
    return sin(M_PI*x)/(M_PI*x);
}

// Window function in Fourier space
double W_NGP(double k)
{
  return sinc( (0.5*k)/GV.KN );
}

double W_CIC(double k)
{
  double sinceval = sinc( (0.5*k)/GV.KN );
  return POW2(sinceval);
}

double W_TSC(double k)
{
  double sinceval = sinc( (0.5*k)/GV.KN );
  return POW3(sinceval);
}

double W_D20_Re(double k)
{
  return gsl_spline_eval(splineRe, (M_PI*k)/GV.KN, acc);
}

double W_D20_Im(double k)
{
  return gsl_spline_eval(splineIm, (M_PI*k)/GV.KN, acc);
}

double zero(double k)
{
  return 0.0;
}


//Sum square window function
double Sum_W2_NGP(double k)
{
  return 1.0;
}

double Sum_W2_CIC(double k)
{
  double sine = sin( (M_PI*0.5*k)/GV.KN );
  return 1.0 - 0.666666666666667*POW2(sine);
}

double Sum_W2_TSC(double k)
{
  double sine = sin( (M_PI*0.5*k)/GV.KN );
  return 1.0 - POW2(sine) + 0.133333333333333*POW4(sine);
}

double Sum_W2_D20(double k)
{
  return 1.0;
}
