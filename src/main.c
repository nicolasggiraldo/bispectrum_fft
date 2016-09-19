#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_spline.h>
#include <fftw3.h>
#include "macros_and_structures.h"
#include "functions.h"



//////////////////////
// GLOBAL VARIABLES //
//////////////////////
/* Global variables for the FFTW routines */
fftw_complex *denConX=NULL; // Density contrast in X-space
fftw_complex *denConK=NULL; // Density contrast in K-space
fftw_plan forwardPlan;
/* Struct Global Variables */
struct globalVariables GV;
/* Global variables for interpolation in Daubechies-type schemes */
gsl_interp_accel *acc=NULL;
gsl_spline *splineRe =NULL;
gsl_spline *splineIm =NULL;
gsl_spline *splineMag2=NULL;
int len_array_D20 = 0;
double *k_D20 =NULL;
double *Re_D20=NULL;
double *Im_D20=NULL;
double *Kmag2_D20=NULL;
FILE *fin_D20=NULL;

typedef double realType;

int main(int argc, char *argv[])
{
  int i,j,k;
  long longi;
  long double Bklongdouble, Ntrilongdouble;

  double (*W_k_Re)(double)=NULL; // Addres to the window function
  double (*W_k_Im)(double)=NULL; // Addres to the window function
  double (*W2_k)(double)  =NULL; // Addres to the window function

  double step=0.5; // Delta_K step for binning
  
  FILE *fout=NULL; // File handler to output
  
  // FFTW complex variables
  fftw_complex *complex=NULL;
  fftw_complex *unity  =NULL;
  realType *real1=NULL;
  realType *real2=NULL;
  realType *real3=NULL;
  realType *I1   =NULL;
  realType *I2   =NULL;
  realType *I3   =NULL;
  fftw_plan realPlan;
  fftw_plan NtriPlan;
  int n[3]; // Number of grids in each axis for FFT estimation
  double *kpos =NULL; // Positions array values according to FFTW k-position convention
  int *indexpos=NULL; // Index values array according to FFTW k-position convention
  double kMag; // Magnitude of the wavevector

  double denConCor[2];
  double WRe, WIm, WMag2;
  double D20_complexki0,D20_complexkj0,D20_complexkk0;
  double D20_complexki1,D20_complexkj1,D20_complexkk1;
  
  
  
  //////////////////////////////
  //* READING PARAMETER FILE *//
  //////////////////////////////
  if(argc<2)
    {
      printf("\n***********************************");
      printf("***********************************\n");
      printf("%s: You must specify the name of the parameter file\n",argv[0]);
      printf("For example: %s pathOfFile/parameterFile.txt\n",argv[0]);
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    }

  // Reading parameter file and verifying there is no error.
  switch( read_parameters(argv[1]) )
    {
    case -1 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: Bad path (or name) to the parameter file.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    case -2 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: Bad settings in the parameter file.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    }
  
  
  
  ///////////////////////////////
  //* READING CELL BINARY FILE *//
  ////////////////////////////////
  /* Read binary cell calcula lee el archivo binario con las
     densidades, despues carga los datos en las variables globales
     y luego calcula la transformada de fourier in-place */
  switch( readBinaryFile() )
    {
    case -1 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: The parameter file could not be allocated.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    case -2 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: FFTW arrays could not be allocated.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    }
  
  
  
  ///////////////////////////////////////
  //* SETTING FFTW3 COORDINATE SYSTEM *//
  ///////////////////////////////////////
  
  /* Position array for storing in the densityContrast */
  kpos = (double *) calloc(GV.NGRID, sizeof(double));
  if(kpos == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: kpos array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  indexpos = (int *) calloc(GV.NGRID, sizeof(int));
  if(indexpos == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: indexpos array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  /* Setting index and positions according to FFTW convention  */
  for( i=0; i<GV.NGRID; i++ ){
    /* REMEMBER kF is (2.0*PI)/L */
    kpos[i]     = (i<GV.NGRID/2) ? GV.KF*i : GV.KF*(i-GV.NGRID);
    indexpos[i] = (i<GV.NGRID/2) ?       i :       (i-GV.NGRID);
  }// for i
  
  
  
  ///////////////////////////////////////////
  //* DEFINING THE WINDOW FUNCTION TO USE *//
  ///////////////////////////////////////////
  if( strcmp(GV.SCHEME, "NGP") == 0 )
    {
      //NGP
      W_k_Re = W_NGP;
      W_k_Im = zero;
      W2_k = Sum_W2_NGP;
    }
  else if( strcmp(GV.SCHEME, "CIC") == 0 )
    {
      //CIC
      W_k_Re = W_CIC; 
      W_k_Im = zero;
      W2_k = Sum_W2_CIC;
    }
  else if( strcmp(GV.SCHEME, "TSC") == 0 )
    {
      // TSC
      W_k_Re = W_TSC; 
      W_k_Im = zero;
      W2_k = Sum_W2_TSC;
    }
  else if( strcmp(GV.SCHEME, "D20") == 0 )
    {
      
      len_array_D20 = 399;
      int err;
      
      k_D20     = (double *) calloc( len_array_D20, sizeof(double));
      Re_D20    = (double *) calloc( len_array_D20, sizeof(double));
      Im_D20    = (double *) calloc( len_array_D20, sizeof(double));
      Kmag2_D20 = (double *) calloc( len_array_D20, sizeof(double));
      
      fin_D20 = fopen("files/D20.txt", "r");
      
      for(i=0; i<len_array_D20; i++)
	{
	  err=fscanf( fin_D20, "%lf %lf %lf %lf",
		      &k_D20[i], &Re_D20[i], &Im_D20[i], &Kmag2_D20[i] );
	}

      fclose(fin_D20);
      
      // GSL interpolation allocation
      acc = gsl_interp_accel_alloc(); // accelerator
      
      // spline interpolation
      splineRe   = gsl_spline_alloc( gsl_interp_cspline, len_array_D20);
      splineIm   = gsl_spline_alloc( gsl_interp_cspline, len_array_D20);
      splineMag2 = gsl_spline_alloc( gsl_interp_cspline, len_array_D20);
      
      // GSL init
      gsl_spline_init(splineRe,   k_D20, Re_D20,    len_array_D20);
      gsl_spline_init(splineIm,   k_D20, Im_D20,    len_array_D20);
      gsl_spline_init(splineMag2, k_D20, Kmag2_D20, len_array_D20);
      
      if(err){}

      // D20
      W_k_Re = zero;
      W_k_Im = zero;
      W2_k = Sum_W2_D20;
    }
  
  
  
  //////////////////////////////
  //* BISPECTRUM CALCULATION *//
  //////////////////////////////
  // Memory allocation
  complex = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  unity   = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  real1   = (realType *) calloc( GV.NGRID3, sizeof(realType));
  real2   = (realType *) calloc( GV.NGRID3, sizeof(realType));
  real3   = (realType *) calloc( GV.NGRID3, sizeof(realType));
  I1      = (realType *) calloc( GV.NGRID3, sizeof(realType));
  I2      = (realType *) calloc( GV.NGRID3, sizeof(realType));
  I3      = (realType *) calloc( GV.NGRID3, sizeof(realType));

  n[X] = n[Y] = n[Z] = GV.NGRID;

  /* Generando el plan  */
  realPlan = fftw_plan_dft(3, n, complex, complex, FFTW_BACKWARD, FFTW_MEASURE);
  NtriPlan = fftw_plan_dft(3, n, unity,   unity,   FFTW_BACKWARD, FFTW_MEASURE);
  
  
  
  /////////////
  // PARA k1 //
  /////////////
  GV.Pk1 = 0.0;
  GV.Nk1 = 0;
  for(i=0; i<GV.NGRID; i++)
    {
      for(j=0; j<GV.NGRID; j++)
	{
	  for(k=0; k<GV.NGRID; k++)
	    {
	      
	      kMag = VECTORMAG(kpos[i],kpos[j],kpos[k]);

	      longi = INDEX(i,j,k);
	      
	      if( (GV.K1-GV.DELTA_K*0.5 < kMag) && (kMag < GV.K1+GV.DELTA_K*0.5) )
		{
		  
		  if( strcmp(GV.SCHEME, "D20") == 0 ) // D20
		    { 
		      W_D20(kpos[i], &D20_complexki0, &D20_complexki1);
		      W_D20(kpos[j], &D20_complexkj0, &D20_complexkj1);
		      W_D20(kpos[k], &D20_complexkk0, &D20_complexkk1);

		      WRe = ( +D20_complexki0 * D20_complexkj0 * D20_complexkk0
			      -D20_complexki0 * D20_complexkj1 * D20_complexkk1
			      -D20_complexki1 * D20_complexkj0 * D20_complexkk1
			      -D20_complexki1 * D20_complexkj1 * D20_complexkk0);
		      
		      WIm = ( +D20_complexki0 * D20_complexkj0 * D20_complexkk1
			      +D20_complexki0 * D20_complexkj1 * D20_complexkk0
			      +D20_complexki1 * D20_complexkj0 * D20_complexkk0
			      -D20_complexki1 * D20_complexkj1 * D20_complexkk1);
		      
		      WMag2 = (WRe*WRe) + (WIm*WIm);
		      
		      denConCor[0] = (denConK[longi][0]*WRe + denConK[longi][1]*WIm)/WMag2;
		      denConCor[1] = (denConK[longi][1]*WRe - denConK[longi][0]*WIm)/WMag2;
		      
		    }
		  else
		    { // NGP, CIC, TSC
		      
		      denConCor[0] = denConK[longi][0]/( W_k_Re(kpos[i]) * W_k_Re(kpos[j]) * W_k_Re(kpos[j]) );
		      denConCor[1] = denConK[longi][1]/( W_k_Re(kpos[i]) * W_k_Re(kpos[j]) * W_k_Re(kpos[j]) );
		      
		    }

		  complex[longi][0] = denConCor[0];
		  complex[longi][1] = denConCor[1];

		  //complex[longi][0] = denConK[longi][0];
		  //complex[longi][1] = denConK[longi][1];

		  unity[longi][0]   = 1.0;
		  unity[longi][1]   = 0.0;

		  GV.Pk1 += COMPLEXMAG(denConK,longi)/( W2_k(kpos[i])*W2_k(kpos[j])*W2_k(kpos[k]) );
		  GV.Nk1++;
		  
		}
	      else
		{
		  
		  complex[longi][0] = 0.0;
		  complex[longi][1] = 0.0;
		  
		  unity[longi][0]   = 0.0;
		  unity[longi][1]   = 0.0;

		}
	    }// for k
	}// for j
    }// for i
  
  /* Do forward FFT */
  // De aqui obtengo real1
  fftw_execute(realPlan);
  fftw_execute(NtriPlan);

  // Realizando copia del resultado...
  for(longi=0; longi<GV.NGRID3; longi++)
    {
      real1[longi] = complex[longi][0];
      I1[longi]    = unity[longi][0];
    }

  // Stimating power spectrum in k1
  GV.Pk1 /= (1.0*GV.Nk1);
  GV.Pk1 *= (GV.SIM_VOL / (1.0 * GV.NGRID3));
  GV.Pk1 *= (       1.0 / (1.0 * GV.NGRID3));
  GV.Pk1 -= GV.SHOT_NOISE;
  GV.Pk1_Error = (GV.Pk1*(GV.DELTA_K/GV.K1))/sqrt(2.0*M_PI);
  
  
  
  /////////////
  // PARA k2 //
  /////////////
  /* Generando el plan  */
  GV.Pk2 = 0.0;
  GV.Nk2 = 0;
  for(i=0; i<GV.NGRID; i++)
    {
      for(j=0; j<GV.NGRID; j++)
	{
	  for(k=0; k<GV.NGRID; k++)
	    {
	      
	      kMag = VECTORMAG(kpos[i],kpos[j],kpos[k]);

	      longi = INDEX(i,j,k);
	      
	      if( (GV.K2-GV.DELTA_K*0.5 < kMag) && (kMag < GV.K2+GV.DELTA_K*0.5) )
		{

		  if( strcmp(GV.SCHEME, "D20") == 0 ) // D20
		    { 
		      W_D20(kpos[i], &D20_complexki0, &D20_complexki1);
		      W_D20(kpos[j], &D20_complexkj0, &D20_complexkj1);
		      W_D20(kpos[k], &D20_complexkk0, &D20_complexkk1);
		      
		      WRe = ( +D20_complexki0 * D20_complexkj0 * D20_complexkk0
			      -D20_complexki0 * D20_complexkj1 * D20_complexkk1
			      -D20_complexki1 * D20_complexkj0 * D20_complexkk1
			      -D20_complexki1 * D20_complexkj1 * D20_complexkk0);
		      
		      WIm = ( +D20_complexki0 * D20_complexkj0 * D20_complexkk1
			      +D20_complexki0 * D20_complexkj1 * D20_complexkk0
			      +D20_complexki1 * D20_complexkj0 * D20_complexkk0
			      -D20_complexki1 * D20_complexkj1 * D20_complexkk1);
		      
		      WMag2 = (WRe*WRe) + (WIm*WIm);
		      
		      denConCor[0] = (denConK[longi][0]*WRe + denConK[longi][1]*WIm)/WMag2;
		      denConCor[1] = (denConK[longi][1]*WRe - denConK[longi][0]*WIm)/WMag2;
		      
		    }
		  else
		    { // NGP, CIC, TSC
		      
		      denConCor[0] = denConK[longi][0]/( W_k_Re(kpos[i]) * W_k_Re(kpos[j]) * W_k_Re(kpos[j]) );
		      denConCor[1] = denConK[longi][1]/( W_k_Re(kpos[i]) * W_k_Re(kpos[j]) * W_k_Re(kpos[j]) );
		      
		    }

		  complex[longi][0] = denConCor[0];
		  complex[longi][1] = denConCor[1];
		  
		  //complex[longi][0] = denConK[longi][0];
		  //complex[longi][1] = denConK[longi][1];

		  unity[longi][0]   = 1.0;
		  unity[longi][1]   = 0.0;

		  GV.Pk2 += COMPLEXMAG(denConK, longi)/( W2_k(kpos[i]) * W2_k(kpos[j]) * W2_k(kpos[k]) );
		  GV.Nk2++;
		  
		}
	      else
		{
		  
		  complex[longi][0] = 0.0;
		  complex[longi][1] = 0.0;
		  
		  unity[longi][0]   = 0.0;
		  unity[longi][1]   = 0.0;
		  
		}
	    }// for k
	}// for j
    }// for i
  
  /* Do forward FFT */
  // De aqui obtengo real2
  fftw_execute(realPlan);
  fftw_execute(NtriPlan);
  
  // Realizando copia del resultado...
  for(longi=0; longi<GV.NGRID3; longi++)
    {
      real2[longi] = complex[longi][0];
      I2[longi]    = unity[longi][0];
    }
  
  // Stimating power spectrum in k2
  GV.Pk2 /= (1.0*GV.Nk2);
  GV.Pk2 *= (GV.SIM_VOL / (1.0 * GV.NGRID3));
  GV.Pk2 *= (       1.0 / (1.0 * GV.NGRID3));
  GV.Pk2 -= GV.SHOT_NOISE;
  GV.Pk2_Error = (GV.Pk2*(GV.DELTA_K/GV.K2))/sqrt(2.0*M_PI);
  
  
  
  /////////////
  // PARA k3 //
  /////////////
  /* Saving data in the outfile */
  printf("Bispectrum calculation\n");
  printf("\n-----------------------------------------------\n");
  printf("Saving data in %s\n", GV.OUTPUT);
  
  fout = fopen(GV.OUTPUT, "w");
  if(fout == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Outfile could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }

  /* Writing header */
  fprintf(fout,"# NGRID          = %d\n",  GV.NGRID);
  fprintf(fout,"# GADGET VERSION = %d\n",  GV.GADGET_VERSION);
  fprintf(fout,"# L              = %g\n",  GV.L);
  fprintf(fout,"# SIM VOL        = %g\n",  GV.SIM_VOL);
  fprintf(fout,"# NP TOT         = %ld\n", GV.NP_TOT);
  fprintf(fout,"# TOTAL MASS     = %lf\n", GV.TOTAL_MASS);
  fprintf(fout,"# RHO MEAN       = %g\n",  GV.RHO_MEAN);
  fprintf(fout,"# VOL_CELL       = %g\n",  GV.VOL_CELL);
  fprintf(fout,"# H              = %g\n",  GV.H);
  fprintf(fout,"# DELTA k        = %g\n",  GV.DELTA_K);
  fprintf(fout,"# kF             = %g\n",  GV.KF);
  fprintf(fout,"# kN             = %lf\n", GV.KN);
  fprintf(fout,"# Shot Noise     = %g\n",  GV.SHOT_NOISE);
  fprintf(fout,"# SCHEME         = %s\n",  GV.SCHEME);
  fprintf(fout,"# OMEGA_M0       = %lf\n", GV.OMEGA_M0);
  fprintf(fout,"# OMEGA_L0       = %lf\n", GV.OMEGA_L0);
  fprintf(fout,"# ZRS            = %lf\n", GV.ZRS);
  fprintf(fout,"# HUBBLEPARAM    = %lf\n", GV.HUBBLEPARAM);
  //fprintf(fout,"\n");
  
  fprintf(fout,"#%19s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n",
	  "k1", "k2", "k3", "P(k1)", "P(k2)", "P(k3)", "B(k1,k2,k3)", "Q(k1,k2,k3)",
	  "Error_P(k1)", "Error_P(k2)", "Error_P(k3)", "Error_B123", "Error_Q123");
  
  /* Lineal binning */
  GV.K3 = GV.DELTA_K*0.5;
  while(GV.K3 <= GV.KN)
    {
      printf("\nk3 = %12.6lf ", GV.K3);

      if( fabs(GV.K1-GV.K2)-(1.5*GV.DELTA_K)<=GV.K3 &&  GV.K3<=(GV.K1+GV.K2)+(1.5*GV.DELTA_K) )
	{
	  printf("skipped!");
	  GV.K3 += (step*GV.DELTA_K);
	  continue;
	}
      
      GV.Pk3 = 0.0;
      GV.Nk3 = 0;
      for(i=0; i<GV.NGRID; i++)
	{
	  for(j=0; j<GV.NGRID; j++)
	    {
	      for(k=0; k<GV.NGRID; k++)
		{
		  
		  kMag = VECTORMAG(kpos[i],kpos[j],kpos[k]);
		  
		  longi = INDEX(i,j,k);
		  
		  if( (GV.K3-GV.DELTA_K*0.5 < kMag) && (kMag < GV.K3+GV.DELTA_K*0.5) )
		    {

		      if( strcmp(GV.SCHEME, "D20") == 0 ) // D20
			{ 
			  W_D20(kpos[i], &D20_complexki0, &D20_complexki1);
			  W_D20(kpos[j], &D20_complexkj0, &D20_complexkj1);
			  W_D20(kpos[k], &D20_complexkk0, &D20_complexkk1);
			  
			  WRe = ( +D20_complexki0 * D20_complexkj0 * D20_complexkk0
				  -D20_complexki0 * D20_complexkj1 * D20_complexkk1
				  -D20_complexki1 * D20_complexkj0 * D20_complexkk1
				  -D20_complexki1 * D20_complexkj1 * D20_complexkk0);
			  
			  WIm = ( +D20_complexki0 * D20_complexkj0 * D20_complexkk1
				  +D20_complexki0 * D20_complexkj1 * D20_complexkk0
				  +D20_complexki1 * D20_complexkj0 * D20_complexkk0
				  -D20_complexki1 * D20_complexkj1 * D20_complexkk1);
			  
			  WMag2 = (WRe*WRe) + (WIm*WIm);
			  
			  denConCor[0] = (denConK[longi][0]*WRe + denConK[longi][1]*WIm)/WMag2;
			  denConCor[1] = (denConK[longi][1]*WRe - denConK[longi][0]*WIm)/WMag2;
			  
			}
		      else
			{ // NGP, CIC, TSC
			  
			  denConCor[0] = denConK[longi][0]/( W_k_Re(kpos[i]) * W_k_Re(kpos[j]) * W_k_Re(kpos[j]) );
			  denConCor[1] = denConK[longi][1]/( W_k_Re(kpos[i]) * W_k_Re(kpos[j]) * W_k_Re(kpos[j]) );
			  
			}

		      complex[longi][0] = denConCor[0];
		      complex[longi][1] = denConCor[1];
		      
		      //complex[longi][0] = denConK[longi][0];
		      //complex[longi][1] = denConK[longi][1];
		      
		      unity[longi][0]   = 1.0;
		      unity[longi][1]   = 0.0;
		      
		      GV.Pk3 += COMPLEXMAG(denConK, longi)/( W2_k(kpos[i]) * W2_k(kpos[j]) * W2_k(kpos[k]) );
		      GV.Nk3++;
		      
		    }
		  else
		    {
		      
		      complex[longi][0] = 0.0;
		      complex[longi][1] = 0.0;
		      
		      unity[longi][0]   = 0.0;
		      unity[longi][1]   = 0.0;
		      
		    }
		}// for k
	    }// for j
	}// for i
      
      /* Do forward FFT */
      // De aqui obtengo real3
      fftw_execute(realPlan);
      fftw_execute(NtriPlan);

      // Realizando copia del resultado...
      for(longi=0; longi<GV.NGRID3; longi++)
	{
	  real3[longi] = complex[longi][0];
	  I3[longi]    = unity[longi][0];
	}

      // Stimating power spectrum in k3
      GV.Pk3 /= (1.0*GV.Nk3);
      GV.Pk3 *= (GV.SIM_VOL / (1.0 * GV.NGRID3));
      GV.Pk3 *= (       1.0 / (1.0 * GV.NGRID3));
      GV.Pk3 -= GV.SHOT_NOISE;
      GV.Pk3_Error = (GV.Pk3*(GV.DELTA_K/GV.K3))/sqrt(2.0*M_PI);
      
      /* Suming values of Bk and Ntri  */
      Bklongdouble = 0.0;
      Ntrilongdouble = 0.0;
      //GV.Bk   = 0.0;
      //GV.Ntri = 0.0;
      for(longi=0; longi<GV.NGRID3; longi++)
	{
	  Bklongdouble  +=(real1[longi]*real2[longi]*real3[longi]);
	  Ntrilongdouble+=(I1[longi]   *I2[longi]   *I3[longi])   ;
	  //GV.Bk   += (real1[longi] * real2[longi] * real3[longi])/(1.0*GV.NGRID3);
	  //GV.Ntri += (I1[longi]    * I2[longi]    * I3[longi])   /(1.0*GV.NGRID3);
	}
      //GV.Bk *= (1.0/GV.NGRID3);
      //GV.Ntri *= (1.0/GV.NGRID3);
      //GV.Ntri  = 8.0*M_PI*M_PI*GV.K1*GV.K2*GV.K3;
      //GV.Ntri *= (GV.DELTA_K/(GV.KF*GV.KF));
      //GV.Ntri *= (GV.DELTA_K/(GV.KF*GV.KF));
      //GV.Ntri *= (GV.DELTA_K/(GV.KF*GV.KF));
      //GV.Bk *= (       1.0/GV.Ntri);

      GV.Bk   = (double)(Bklongdouble/Ntrilongdouble);
      GV.Ntri = (double)(Ntrilongdouble/(1.0*GV.NGRID3));
      
      GV.Bk *= (GV.SIM_VOL/GV.NGRID3);
      GV.Bk *= (GV.SIM_VOL/GV.NGRID3);
      GV.Bk *= (       1.0/GV.NGRID3);
      

      if(GV.Bk<0 || GV.Ntri<1.0)
	{
	  printf("skipped!");
	  GV.K3 += (step*GV.DELTA_K);
	  continue;
	}

      // Stimating shot noise for bispectrum
      GV.Bk_shotnoise = (GV.Pk1+GV.Pk2+GV.Pk3)*GV.SHOT_NOISE + POW2(GV.SHOT_NOISE);
      
      GV.Bk -= GV.Bk_shotnoise;

      GV.Bk_Error  = sqrt(M_PI/(GV.K1*GV.K2*GV.K3*POW3(GV.S_KF)));
      GV.Bk_Error *= sqrt(GV.Pk1*GV.Pk2*GV.Pk3);
      
      GV.Qk = (GV.Bk) / (GV.Pk1*GV.Pk2+GV.Pk2*GV.Pk3+GV.Pk1*GV.Pk3);

      GV.Qk_Error  = GV.Bk_Error;
      GV.Qk_Error -= GV.Qk*( (GV.Pk2+GV.Pk3)*GV.Pk1_Error +
			     (GV.Pk3+GV.Pk1)*GV.Pk2_Error +
			     (GV.Pk1+GV.Pk2)*GV.Pk3_Error);
      GV.Qk_Error /= ( GV.Pk1*GV.Pk2 + GV.Pk2*GV.Pk3 + GV.Pk3*GV.Pk1 );
      
      fprintf(fout,"%20lf %20lf %20lf %20e %20e %20e %20e %20e %20e %20e %20e %20e %20e\n",
	      GV.K1, GV.K2, GV.K3, GV.Pk1, GV.Pk2, GV.Pk3, GV.Bk, GV.Qk,
	      GV.Pk1_Error, GV.Pk2_Error, GV.Pk3_Error,
	      GV.Bk_Error, GV.Qk_Error);
      
      GV.K3 += (step*GV.DELTA_K);
      //GV.K3 += (0.5 * GV.DELTA_K);
    
    }// while
  
  
  
  //////////////////////
  //* FREEING MEMORY *//
  //////////////////////
  printf("\nMemory free\n");
  if( strcmp(GV.SCHEME, "D20") == 0 )
    {
      gsl_spline_free( splineRe );
      gsl_spline_free( splineIm );
      gsl_spline_free( splineMag2 );
      gsl_interp_accel_free( acc );
    
      free(k_D20);
      free(Re_D20);
      free(Im_D20);
      free(Kmag2_D20);
    }
  fclose(fout);
  free(GV.FILE_NAME);
  free(GV.OUTPUT);
  fftw_free(denConK);
  fftw_free(complex);
  fftw_free(unity);
  free(real1);
  free(real2);
  free(real3);
  free(I1);
  free(I2);
  free(I3);
  fftw_destroy_plan(realPlan);
  fftw_destroy_plan(NtriPlan);

  if(W_k_Im==NULL || W_k_Re==NULL){}		 
   
  return 0;
}
