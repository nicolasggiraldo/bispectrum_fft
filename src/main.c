#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_spline.h>
#include <fftw3.h>
#include "macros_and_structures.h"
#include "functions.h"

int readBinaryFile();

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



int main(int argc, char *argv[])
{
  int i,j,k;
  long longi;

  double (*W_k_Re)(double)=NULL; // Addres to the window function
  double (*W_k_Im)(double)=NULL; // Addres to the window function
  double (*W2_k)(double)  =NULL; // Addres to the window function
  
  FILE *fout=NULL; // File handler to output
  
  // FFTW complex variables
  fftw_complex *complex=NULL;
  fftw_complex *real1  =NULL;
  fftw_complex *real2  =NULL;
  fftw_complex *real3  =NULL;
  fftw_complex *unity  =NULL;
  fftw_complex *I1     =NULL;
  fftw_complex *I2     =NULL;
  fftw_complex *I3     =NULL;
  fftw_plan realPlan;
  fftw_plan NtriPlan;
  int n[3]; // Number of grids in each axis for FFT estimation
  double *kpos =NULL; // Positions array values according to FFTW k-position convention
  int *indexpos=NULL; // Index values array according to FFTW k-position convention
  double kMag; // Magnitude of the wavevector

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
      W_k_Re = W_D20_Re;
      W_k_Im = W_D20_Im;
      W2_k = Sum_W2_D20;
    }
  
  
  
  //////////////////////////////
  //* BISPECTRUM CALCULATION *//
  //////////////////////////////
  // Memory allocation
  //complex = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  real1   = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  real2   = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  real3   = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  //unity   = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  I1      = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  I2      = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  I3      = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);

  n[X] = n[Y] = n[Z] = GV.NGRID;
  
  
  
  /////////////
  // PARA k1 //
  /////////////
  /* Generando el plan  */
  complex = real1;
  realPlan = fftw_plan_dft(3, n, complex, real1, FFTW_BACKWARD, FFTW_ESTIMATE);
  I1 = unity;
  NtriPlan = fftw_plan_dft(3, n, unity,   I1,    FFTW_BACKWARD, FFTW_ESTIMATE);

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
		  
		  complex[longi][0] = denConK[longi][0];
		  complex[longi][1] = denConK[longi][1];

		  unity[longi][0]   = 1.0;
		  unity[longi][1]   = 0.0;

		  GV.Pk1 += COMPLEXMAG(denConK, longi)/( W2_k(kpos[i]) * W2_k(kpos[j]) * W2_k(kpos[k]) );
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

  // Destruyendo el plan
  fftw_destroy_plan(realPlan);
  fftw_destroy_plan(NtriPlan);

  // Stimating power spectrum in k1
  GV.Pk1 /= GV.Nk1;
  GV.Pk1 *= (GV.SIM_VOL / (1.0 * GV.NGRID3 * GV.NGRID3));
  GV.Pk1 -= GV.SHOT_NOISE;
  GV.Pk1_Error = (GV.Pk1*(GV.DELTA_K/GV.K1))/sqrt(2.0*M_PI);
  
  
  
  /////////////
  // PARA k2 //
  /////////////
  /* Generando el plan  */
  complex = real2;
  realPlan = fftw_plan_dft(3, n, complex, real2, FFTW_BACKWARD, FFTW_ESTIMATE);
  I2 = unity;
  NtriPlan = fftw_plan_dft(3, n, unity,   I2,    FFTW_BACKWARD, FFTW_ESTIMATE);

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
		  
		  complex[longi][0] = denConK[longi][0];
		  complex[longi][1] = denConK[longi][1];

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

  // Destruyendo el plan
  fftw_destroy_plan(realPlan);
  fftw_destroy_plan(NtriPlan);

  // Stimating power spectrum in k1
  GV.Pk2 /= GV.Nk2;
  GV.Pk2 *= (GV.SIM_VOL / (1.0 * GV.NGRID3 * GV.NGRID3));
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
  fprintf(fout,"# L              = %lf\n", GV.L);
  fprintf(fout,"# SIM VOL        = %lf\n", GV.SIM_VOL);
  fprintf(fout,"# NP TOT         = %ld\n", GV.NP_TOT);
  fprintf(fout,"# TOTAL MASS     = %lf\n", GV.TOTAL_MASS);
  fprintf(fout,"# RHO MEAN       = %lf\n", GV.RHO_MEAN);
  fprintf(fout,"# VOL_CELL       = %lf\n", GV.VOL_CELL);
  fprintf(fout,"# H              = %lf\n", GV.H);
  fprintf(fout,"# DELTA k        = %lf\n", GV.DELTA_K);
  fprintf(fout,"# kF             = %lf\n", GV.KF);
  fprintf(fout,"# kN             = %lf\n", GV.KN);
  fprintf(fout,"# Shot Noise     = %lf\n", GV.SHOT_NOISE);
  fprintf(fout,"# SCHEME         = %s\n",  GV.SCHEME);
  fprintf(fout,"# OMEGA_M0       = %lf\n", GV.OMEGA_M0);
  fprintf(fout,"# OMEGA_L0       = %lf\n", GV.OMEGA_L0);
  fprintf(fout,"# ZRS            = %lf\n", GV.ZRS);
  fprintf(fout,"# HUBBLEPARAM    = %lf\n", GV.HUBBLEPARAM);
  fprintf(fout,"\n");
  
  fprintf(fout,"#%19s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n",
	  "k1", "k2", "k3", "P(k1)", "P(k2)", "P(k3)", "B(k1,k2,k3)", "Q(k1,k2,k3)",
	  "Error_P(k1)", "Error_P(k2)", "Error_P(k3)", "Error_B123", "Error_Q123");
  
  /* Generando el plan  */
  complex = real3;
  realPlan = fftw_plan_dft(3, n, complex, real3, FFTW_BACKWARD, FFTW_MEASURE);
  I3 = unity;
  NtriPlan = fftw_plan_dft(3, n, unity,   I3,    FFTW_BACKWARD, FFTW_MEASURE);
  
  
  /* Lineal binning */
  GV.K3 = GV.DELTA_K*0.5;
  while(GV.K3 <= GV.KN)
    {
      printf("k3 = %g\n", GV.K3);
      
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
		      
		      complex[longi][0] = denConK[longi][0];
		      complex[longi][1] = denConK[longi][1];
		      
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
      
      // Stimating power spectrum in k3
      GV.Pk3 /= GV.Nk3;
      GV.Pk3 *= (GV.SIM_VOL / (1.0 * GV.NGRID3 * GV.NGRID3));
      GV.Pk3 -= GV.SHOT_NOISE;
      GV.Pk3_Error = (GV.Pk3*(GV.DELTA_K/GV.K3))/sqrt(2.0*M_PI);
      
      /* Suming values of Bk and Ntri  */
      GV.Bk   = 0.0;
      GV.Ntri = 0.0;
      for(longi=0; longi<GV.NGRID3; longi++)
	{
	  GV.Bk   += (real1[longi][0] * real2[longi][0] * real3[longi][0]);
	  GV.Ntri += (I1[longi][0]    * I2[longi][0]    * I3[longi][0]);
	}
      
      GV.Bk *= (GV.SIM_VOL/GV.NGRID3);
      GV.Bk *= (GV.SIM_VOL/GV.NGRID3);
      GV.Bk *= (      1.0/GV.NGRID3);
      GV.Bk *= (      1.0/GV.Ntri);
      
      // Stimating shot noise for bispectrum
      GV.Bk_shotnoise = (GV.Pk1+GV.Pk2+GV.Pk3)*GV.SHOT_NOISE + POW2(GV.SHOT_NOISE);
      
      GV.Bk -= GV.Bk_shotnoise;

      GV.Bk_Error  = sqrt(M_PI/(GV.K1*GV.K2*GV.K3*POW3(GV.S_KF)));
      GV.Bk_Error *= sqrt(GV.Pk1*GV.Pk2*GV.Pk3);
      
      GV.Qk = GV.Bk / (GV.Pk1*GV.Pk2+GV.Pk2*GV.Pk3+GV.Pk1*GV.Pk3);

      GV.Qk_Error  = GV.Bk_Error;
      GV.Qk_Error -= GV.Qk*( (GV.Pk2+GV.Pk3)*GV.Pk1_Error +
			     (GV.Pk3+GV.Pk1)*GV.Pk2_Error +
			     (GV.Pk1+GV.Pk2)*GV.Pk3_Error);
      GV.Qk_Error /= ( GV.Pk1*GV.Pk2 + GV.Pk2*GV.Pk3 + GV.Pk3*GV.Pk1 );
      
      fprintf(fout,"%20lf %20lf %20lf %20e %20e %20e %20e %20e %20e %20e %20e %20e %20e\n",
	      GV.K1, GV.K2, GV.K3, GV.Pk1, GV.Pk2, GV.Pk3, GV.Bk, GV.Qk,
	      GV.Pk1_Error, GV.Pk2_Error, GV.Pk3_Error,
	      GV.Bk_Error, GV.Qk_Error);
      
      GV.K3 += (GV.DELTA_K);
      //GV.K3 += (0.5 * GV.DELTA_K);
    
    }// while
  

  /* Freeing up memory allocation */
  fclose(fout);
  free(GV.FILE_NAME);
  free(GV.OUTPUT);
  fftw_free(denConK);
  fftw_free(complex);
  fftw_free(real1);
  fftw_free(real2);
  fftw_free(real3);
  fftw_free(unity);
  fftw_free(I1);
  fftw_free(I2);
  fftw_free(I3);
  fftw_destroy_plan(realPlan);
  fftw_destroy_plan(NtriPlan);
   
  return 0;
}



/*
 * Function:  readBinaryFile
 * --------------------                  
 * Reads a binary file with the information of the cell and store 
 * the information in the data structure variable *part* it also 
 * returns the total number of particles.
 * 
 *  There are no arguments in the routiene.              
 *
 *  returns: Integer value.               
 *            0 --> There is no error. 
 *           -1 --> There is an error loading file
 *           -2 --> Structure cell could not be allocated. 
 */
int readBinaryFile()
{
  
  FILE *fdata;
  long id_cell;
  double density_Contrast;
  int n[3]; // Number of grids in each axis for FFT estimation
  size_t err;
  
  printf("\n-----------------------------------------------\n");
  printf("Reading file:   %s\n", GV.FILE_NAME);
  
  fdata = fopen(GV.FILE_NAME,"rb");
  if(fdata == NULL)
    {
      printf("File %s cannot be open\n", GV.FILE_NAME);
      return -1;
    }

  /* Getting cosmological parameters of the simulation */
  err = fread(&GV.OMEGA_M0,    sizeof(double), 1, fdata);
  err = fread(&GV.OMEGA_L0,    sizeof(double), 1, fdata);
  err = fread(&GV.ZRS,         sizeof(double), 1, fdata);
  err = fread(&GV.HUBBLEPARAM, sizeof(double), 1, fdata);

  /* Getting simulation parameters */
  err = fread(&GV.NGRID,          sizeof(int),      1, fdata);
  err = fread(&GV.GADGET_VERSION, sizeof(int),      1, fdata);
  err = fread(&GV.L,              sizeof(double),   1, fdata);
  err = fread(&GV.NP_TOT,         sizeof(long int), 1, fdata);
  err = fread(&GV.TOTAL_MASS,     sizeof(double),   1, fdata);
  err = fread(&GV.RHO_MEAN,       sizeof(double),   1, fdata);
  err = fread(&GV.VOL_CELL,       sizeof(double),   1, fdata);
  err = fread(&GV.H,              sizeof(double),   1, fdata);
  err = fread(&(GV.SCHEME[0]),    sizeof(char),     1, fdata);
  err = fread(&(GV.SCHEME[1]),    sizeof(char),     1, fdata);
  err = fread(&(GV.SCHEME[2]),    sizeof(char),     1, fdata);
  GV.SCHEME[3] = '\0';

  GV.NGRID2     = POW2(GV.NGRID);
  GV.NGRID3     = POW3(GV.NGRID);
  GV.SIM_VOL    = POW3(GV.L);
  GV.KF         = (2.0*M_PI) / GV.L;
  GV.DELTA_K    = GV.S_KF * GV.KF;
  GV.SHOT_NOISE = GV.VOL_CELL / GV.NP_TOT;
  GV.KN         = M_PI / GV.H;

  printf("\n-----------------------------------------------\n");
  printf("The original snapshot has a total of %ld particles\n", GV.NP_TOT);
  printf("----------------------------------------\n");
  printf(" * Redshift...     %20.8lf\n", GV.ZRS);
  printf(" * Omega0...       %20.8lf\n", GV.OMEGA_M0);
  printf(" * OmageLa...      %20.8lf\n", GV.OMEGA_L0);
  printf(" * Hubbleparam...  %20.8lf\n", GV.HUBBLEPARAM);
  printf("----------------------------------------\n");
  printf(" * Boxsize...      %20.8lf\n", GV.L);
  printf(" * Ngrid...        %20d\n",    GV.NGRID);
  printf(" * SimMass...      %20.8e\n",  GV.TOTAL_MASS);
  printf(" * Scheme...       %20s\n",    GV.SCHEME);
  printf("----------------------------------------\n");
  printf(" * kF...           %20.8lf\n", GV.KF);
  printf(" * kN...           %20.8lf\n", GV.KN);
  printf(" * DELTA_k...      %20.8lf\n", GV.DELTA_K);
  printf(" * PSshotNoise...  %20.8e\n",  GV.SHOT_NOISE);

  

  /////////////////////////////////
  //* FFTW VARIABLES ALLOCATION *//
  /////////////////////////////////
  denConK = fftw_malloc( sizeof(fftw_complex) * GV.NGRID3 );
  if(denConK == NULL)
    {
      printf("FFTW arrays could not be allocated\n");
      return -2;
    }//if
  
  denConX = fftw_malloc( sizeof(fftw_complex) * GV.NGRID3 );
  if(denConX == NULL)
    {
      printf("FFTW arrays could not be allocated\n");
      return -2;
    }//if
  
  
  
  /////////////////////////
  //* ESTABLISHING PLAN *//
  /////////////////////////
  // Number of grids in each axis
  n[X] = n[Y] = n[Z] = GV.NGRID;
  
  // FFTW plan for a 3D Fourier transform
  /* FFTW_MEASURE: find the optimal plan by actually computing several FFTs 
     FFTW_ESTIMATE: do not run any FFT and provide a "reasonable" plan
     FFTW_OUT_OF_PLACE: a plan assumes that the in and out are distinct 
     FFTW_IN_PLACE: a plan assumes that the in and out are same */
  forwardPlan = fftw_plan_dft(3, n, denConX, denConK, FFTW_FORWARD, FFTW_ESTIMATE);
  
  
  
  ////////////////////////////////
  //* GETTING DENSITY CONTRAST *//
  ////////////////////////////////
  for(id_cell=0; id_cell<GV.NGRID3; id_cell++)
    {
      err = fread(&density_Contrast, sizeof(double), 1, fdata);
      printf("%lf\n", density_Contrast);
      denConX[id_cell][0] = density_Contrast;
      denConX[id_cell][1] = 0.0;
    }//for id_cell
  
  
  
  ////////////////////////////////////
  //* TAKING THE FOURIER TRANSFORM *//
  ////////////////////////////////////
  fftw_execute(forwardPlan);
  printf("\n-----------------------------------------------\n");
  printf("Fourier Transform succes.\n");
  fftw_free(denConX);
  fftw_destroy_plan(forwardPlan);
    
  
  
  if(err){};
  
  fclose(fdata);
  
  return 0;
}
