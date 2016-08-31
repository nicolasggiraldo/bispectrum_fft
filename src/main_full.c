#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

/*
  FFTW_MEASURE: find the optimal plan by actually computing several FFTs 
  FFTW_ESTIMATE: do not run any FFT and provide a "reasonable" plan
  FFTW_OUT_OF_PLACE: a plan assumes that the in and out are distinct 
  FFTW_IN_PLACE: a plan assumes that the in and out are same
*/

//#define DK 0.02
#define TOL 0.01
#define ZERO 1e-7
#define INDEX(i,j,k) (k)+GV.NGRID*((j)+GV.NGRID*(i)) /* Index preprocessor 
                                                        for the C-Order 
						     */
#define COMPLEXMAG(A,i) ( (A[i][0] * A[i][0]) + (A[i][1] * A[i][1]) )
#define VECTORMAG(x,y,z) sqrt( (x * x) + (y * y) + (z * z) )

// Arguments of B(k1, k2, k3) which are going to forma a tringle
double k1, k2, k3;

#include "structures.h"
#include "functions.h"
#include "readWrite.h"

int main(int argc, char *argv[])
{
  
  long i, j, k;
  int n[3];                   // Number of grids in each dimension.
  // FFTW variables
  fftw_complex *denConX=NULL; // Density contrast in X-space
  fftw_complex *denConK=NULL; // Density contrast in K-space
  fftw_plan forwardPlan;

  // For the second method
  fftw_complex *complex=NULL; // Complex
  double *real1=NULL, *real2=NULL, *real3=NULL;
  fftw_plan realPlan;

  // For the calculation of Ntri
  double Ntri;
  fftw_complex *unity=NULL; // Complex
  //fftw_complex *I1=NULL, *I2=NULL, *I3=NULL;
  double *I1=NULL, *I2=NULL, *I3=NULL;
  fftw_plan NtriPlan;
  
  double kx, ky, kz, *kMag=NULL;
  double one_over_NGRID;
  size_t *order=NULL;
  double Bk, Qk; // Bispectrum and reduced bispectrum
  double Bk_shotnoise; 
  double Pk1, Pk2, Pk3;
  // Pointer to outfile
  FILE *fout=NULL;
  

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












  
  read_parameters(argv[1]);
  /* Calculating the number of cells from GV.NGRID */
  GV.NGRID3 = GV.NGRID * GV.NGRID * GV.NGRID;
  GV.dx     = GV.L / GV.NGRID;
  GV.SimVol = GV.L * GV.L * GV.L;
  GV.kF     = (2.0*M_PI)/GV.L;
  GV.kN     = M_PI/GV.dx;
  //GV.deltaK = deltaK();
  GV.deltaK = 3.0 * GV.kF;
  GV.shotNoise = GV.SimVol/GV.NpTot;
  //GV.shotNoise = pow(GV.dx,3);
  one_over_NGRID = 1.0/GV.NGRID;

  printf("kF = %g\n",GV.kF);
  printf("kN = %g\n",GV.kN);
  //GV.kN = round(0.9*GV.kN*100.0)/100.0;
  //printf("New kN = %g\n",GV.kN);
  printf("deltaK = %g\n",GV.deltaK);

  
  ////////////////////////////////
  //* EXECUTING FFTW3 ROUTINES *//
  ////////////////////////////////
  /* Allocate memory for the input and output fftw arrays,
     initialize input array to (1.,0.) */
  denConX = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  denConK = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  /* Routine for reading cell information  */
  readCell(denConX);
  /* Number of grids in each axis */
  n[0] = n[1] = n[2] = GV.NGRID;
  /* FFTW plan for a 3D Fourier transform */
  forwardPlan = fftw_plan_dft(3, n, denConX, denConK,
			      FFTW_FORWARD, FFTW_ESTIMATE);
  
  /* Do forward FFT */
  fftw_execute(forwardPlan);
  printf("Fourier Transform succes!!\n");

  // Freeing up denConX
  fftw_free(denConX);
  
  /* Magnitud of the k vector */
  kMag = (double *)calloc(GV.NGRID3, sizeof(double));

  /* REMEMBER kF is (2.0*PI)/L; */
  for(i=0L; i<GV.NGRID; i++){
    
    // Momentum coordinate in the X-Axis
    kx = (i<GV.NGRID*0.5) ? GV.kF*i : GV.kF*(i-GV.NGRID);
    
    for(j=0L; j<GV.NGRID; j++){
      
      // Momentum coordinate in the Y-Axis
      ky = (j<GV.NGRID*0.5) ? GV.kF*j : GV.kF*(j-GV.NGRID);
      
      for(k=0L; k<GV.NGRID; k++){
	
	// Momentum coordinate in the Z-Axis
	kz = (k<GV.NGRID*0.5) ? GV.kF*k : GV.kF*(k-GV.NGRID);

	// Distance from kX=0, kY=0, kZ=0.
	kMag[INDEX(i,j,k)] = VECTORMAG(kx, ky, kz);
      }// for k
    }// for j
  }// for i
    
  printf("Fourier space positions assigned!!\n");
  
  //////////////////////////////
  //* BISPECTRUM CALCULATION *//
  //////////////////////////////
  
  // Memory allocation
  complex = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  real1   = (double *)calloc(GV.NGRID3, sizeof(double));
  real2   = (double *)calloc(GV.NGRID3, sizeof(double));
  real3   = (double *)calloc(GV.NGRID3, sizeof(double));
  unity   = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  //I1      = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  //I2      = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  //I3      = fftw_malloc(sizeof(fftw_complex)*GV.NGRID3);
  I1      = (double *)calloc(GV.NGRID3, sizeof(double));
  I2      = (double *)calloc(GV.NGRID3, sizeof(double));
  I3      = (double *)calloc(GV.NGRID3, sizeof(double));
  
  /////////////
  // PARA k1 //
  /////////////
  
  /* Generando el plan  */
  realPlan = fftw_plan_dft_c2r(3, n, complex, real1, FFTW_ESTIMATE);
  NtriPlan = fftw_plan_dft_c2r(3, n, unity,   I1,    FFTW_ESTIMATE);
  
  for(i=0L; i<GV.NGRID3; i++){
    if( k1-GV.deltaK*0.5 < kMag[i] && kMag[i] < k1+GV.deltaK*0.5 ){

      complex[i][0] = denConK[i][0];
      complex[i][1] = denConK[i][1];

      unity[i][0]   = 1.0;
      unity[i][1]   = 0.0;
      
    }else{

      complex[i][0] = 0.0;
      complex[i][1] = 0.0;

      unity[i][0]   = 0.0;
      unity[i][1]   = 0.0;
      
    }
  }// for i
  
  /* Do forward FFT */
  // De aqui obtengo real1
  fftw_execute(realPlan);
  fftw_execute(NtriPlan);

  // Destruyendo el plan
  fftw_destroy_plan(realPlan);
  fftw_destroy_plan(NtriPlan);

  // Stimating power spectrum in k1
  /* Suming values of Pk  */
  Pk1 = 0.0;
  for(i=0L;i<GV.NGRID3;i++)
    Pk1 += (real1[i]*real1[i]);

  Pk1 *= ((GV.dx * GV.dx * GV.dx)/(1.0*GV.NGRID3))/(1.0*GV.NGRID3);
  Pk1 *= (GV.kF*GV.kF*GV.kF)/(4.0*M_PI*k1*k1*GV.deltaK);
  Pk1 -= GV.shotNoise;

  
  /////////////
  // PARA k2 //
  /////////////

  /* Generando el plan  */
  realPlan = fftw_plan_dft_c2r(3, n, complex, real2, FFTW_ESTIMATE);
  NtriPlan = fftw_plan_dft_c2r(3, n,   unity,    I2, FFTW_ESTIMATE);

  for(i=0L; i<GV.NGRID3; i++){
    if( k2-GV.deltaK*0.5 < kMag[i] && kMag[i] < k2+GV.deltaK*0.5 ){
      
      complex[i][0] = denConK[i][0];
      complex[i][1] = denConK[i][1];

      unity[i][0]   = 1.0;
      unity[i][1]   = 0.0;
      
    }else{
      
      complex[i][0] = 0.0;
      complex[i][1] = 0.0;

      unity[i][0]   = 0.0;
      unity[i][1]   = 0.0;
      
    }
  }// for i       

  /* Do forward FFT */
  // De aqui obtengo real2
  fftw_execute(realPlan);
  fftw_execute(NtriPlan);

  // Destruyendo el plan
  fftw_destroy_plan(realPlan);
  fftw_destroy_plan(NtriPlan);

  // Stimating power spectrum in k2
  /* Suming values of Pk  */
  Pk2 = 0.0;
  for(i=0L;i<GV.NGRID3;i++)
    Pk2 += (real2[i]*real2[i]);

  Pk2 *= ((GV.dx * GV.dx * GV.dx)/(1.0*GV.NGRID3))/(1.0*GV.NGRID3);
  Pk2 *= (GV.kF*GV.kF*GV.kF)/(4.0*M_PI*k2*k2*GV.deltaK);
  Pk2 -= GV.shotNoise;


  /////////////
  // PARA k3 //
  /////////////

  printf("Bispectrum calculation\n");
  printf("Saving data in %s\n", GV.OUTPUT);
  
  /* Generando el plan  */
  realPlan = fftw_plan_dft_c2r(3, n, complex, real3, FFTW_ESTIMATE);
  NtriPlan = fftw_plan_dft_c2r(3, n,   unity,    I3, FFTW_ESTIMATE);

  /* Numero de bins para los que calcular el power spectrum */
  //Nbins = ceil( GV.kN / GV.deltaK );
  //printf("Nbins=%d\n",Nbins);

  /* Generando salida del programa  */
  fout = fopen(GV.OUTPUT, "w");
  fprintf(fout,"#%19s %20s %20s %20s %20s\n",
	  "k1", "k2", "k3", "B(k1,k2,k3)", "Q(k1,k2,k3)");

  /* Lineal binning */
  k3 = GV.deltaK*0.5;
  while(k3 <= GV.kN){
    //printf("l = %d\n",l);
    printf("k3 = %lf\n", k3);
    //k3 = (l + 0.5) * GV.deltaK;
    for(i=0L; i<GV.NGRID3; i++){
      if( k3-GV.deltaK*0.5 < kMag[i] && kMag[i] < k3+GV.deltaK*0.5 ){
	
	complex[i][0] = denConK[i][0];
	complex[i][1] = denConK[i][1];

	unity[i][0]   = 1.0;
	unity[i][1]   = 0.0;
	
      }else{
	
	complex[i][0] = 0.0;
	complex[i][1] = 0.0;

	unity[i][0]   = 0.0;
	unity[i][1]   = 0.0;
	
      }
    }// for i

    /* Do forward FFT */
    // De aqui obtengo real3
    fftw_execute(realPlan);
    fftw_execute(NtriPlan);

    /* Suming values of Bk and Ntri  */
    Bk = 0.0;
    Ntri = 0.0;
    
    for(i=0L; i<GV.NGRID3; i++){
      Bk   += (real1[i] * real2[i] * real3[i]);
      Ntri += (I1[i]    * I2[i]    * I3[i]);
    }
    
    Ntri *= (1.0/GV.NGRID3);

    Bk *= GV.SimVol*GV.SimVol;
    Bk *= pow(one_over_NGRID, 12);
    //Bk *= (GV.kF*GV.kF*GV.kF*GV.kF*GV.kF*GV.kF)/(8.0*M_PI*M_PI*k1*k2*k3*GV.deltaK*GV.deltaK*GV.deltaK);
    Bk *= (1.0/Ntri);


    // Stimating power spectrum in k3
    /* Suming values of Pk  */
    Pk3 = 0.0;
    for(i=0L;i<GV.NGRID3;i++)
      Pk3 += (real3[i]*real3[i]);
    
    Pk3 *= ((GV.dx * GV.dx * GV.dx)/(1.0*GV.NGRID3))/(1.0*GV.NGRID3);
    Pk3 *= (GV.kF*GV.kF*GV.kF)/(4.0*M_PI*k3*k3*GV.deltaK);
    Pk3 -= GV.shotNoise;

    // Stimating shot noise for bispectrum
    Bk_shotnoise = GV.shotNoise*(Pk1+Pk2+Pk3) + (GV.shotNoise*GV.shotNoise);

    Bk -= Bk_shotnoise;

    Qk = Bk / ( Pk1*Pk2 + Pk2*Pk3 + Pk1*Pk3 );

    fprintf(fout,"%20lf %20lf %20lf %20e %20e\n",k1, k2, k3, Bk, Qk);

    k3 +=0.5 * GV.deltaK;
    //k3 += 0.01;

  }// while

  fclose(fout);
  
  printf("Suerte parce!\n");
  
  /* Freeing up memory allocation */
  free(order);
  free(real1);
  free(real2);
  free(real3);
  free(kMag);
  fftw_free(denConK);
  fftw_free(complex);
  //fftw_free(I1);
  //fftw_free(I2);
  //fftw_free(I3);
  free(I1);
  free(I2);
  free(I3);
  fftw_destroy_plan(forwardPlan);
  fftw_destroy_plan(realPlan);
  fftw_destroy_plan(NtriPlan);
   
  return 0;
}
