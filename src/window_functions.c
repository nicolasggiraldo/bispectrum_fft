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


double h[] = { 2.667005790055555358661744877130858277192498290851289932779975e-02,
	       1.881768000776914890208929736790939942702546758640393484348595e-01,
	       5.272011889317255864817448279595081924981402680840223445318549e-01,
	       6.884590394536035657418717825492358539771364042407339537279681e-01,
	       2.811723436605774607487269984455892876243888859026150413831543e-01,
	       -2.498464243273153794161018979207791000564669737132073715013121e-01,
	       -1.959462743773770435042992543190981318766776476382778474396781e-01,
	       1.273693403357932600826772332014009770786177480422245995563097e-01,
	       9.305736460357235116035228983545273226942917998946925868063974e-02,
	       -7.139414716639708714533609307605064767292611983702150917523756e-02,
	       -2.945753682187581285828323760141839199388200516064948779769654e-02,
	       3.321267405934100173976365318215912897978337413267096043323351e-02,
	       3.606553566956169655423291417133403299517350518618994762730612e-03,
	       -1.073317548333057504431811410651364448111548781143923213370333e-02,
	       1.395351747052901165789318447957707567660542855688552426721117e-03,
	       1.992405295185056117158742242640643211762555365514105280067936e-03,
	       -6.858566949597116265613709819265714196625043336786920516211903e-04,
	       -1.164668551292854509514809710258991891527461854347597362819235e-04,
	       9.358867032006959133405013034222854399688456215297276443521873e-05,
	       -1.326420289452124481243667531226683305749240960605829756400674e-05 };
  

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
  return gsl_spline_eval(splineRe, k/GV.KN, acc);
}

double W_D20_Im(double k)
{
  return gsl_spline_eval(splineIm, k/GV.KN, acc);
}

void Hfun(double w, double *real, double *imag)
{
  int k;

  *real=0.0;
  *imag=0.0;
  
  for(k=0; k<20; k++)
    {    
      *real +=  h[k]*cos(w*(k*1.0));
      *imag += -h[k]*sin(w*(k*1.0));
    }

  *real *= M_SQRT1_2;
  *imag *= M_SQRT1_2;
}

void W_D20(double k, double *WRe, double *WIm)
{
  double arg = (M_PI*k)/GV.KN;
  double prod0[] = {1.0, 0.0};
  double Heval[2], prod1[2];
  
  while(1)
    {
      arg *= 0.5;
      Hfun( arg, &(Heval[0]), &(Heval[1]) );
      prod1[0] = prod0[0]*Heval[0]-prod0[1]*Heval[1];
      prod1[1] = prod0[0]*Heval[1]+prod0[1]*Heval[0];
      
      if( fabs(prod0[0]-prod1[0])<TOL && fabs(prod0[1]-prod1[1])<TOL )
	break;
      
      prod0[0] = prod1[0];
      prod0[1] = prod1[1];
    }				
    
    *WRe=prod1[0];
    *WIm=prod1[1];		
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
