/*
* Author:Pruthvi Suryadevara
* Email: pruthvi.suryadevara@tifr.res.in
* Description: Fourier Transform using GSL
* Compile using gcc P_3.c -lgsl -lgslcblas -lm -o P_3.out
* Note: Instead of multipiying the solution by phase factor we are defining y that is already shifted
* which mwke no difference as phase in fourier plane is shift in x plane
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z, i) ((z)[2 * (i)])
#define IMAG(z, i) ((z)[2 * (i) + 1])


int main()
{
  gsl_fft_complex_wavetable * wavetable;
  gsl_fft_complex_workspace * workspace;
  double *y, *ky, *kys;
  int N = 1024;
  double x, k, c_fac; 
  double range[] = {-100, 100};
  double dx = (range[1] - range[0]) / N;
  double dk = 1 / (N * dx);
  y = (double*) malloc(sizeof(double) * 2 * N);
  ky = (double*) malloc(sizeof(double) * 2 * N);
  kys = (double*) malloc(sizeof(double) * 2 * N);

  REAL(y,0) = 1;
  REAL(ky,0) = 1;
  for (int i=1; i<N/2; i++)   // Defining shifted y 
    {
      x = (i * dx);
      REAL(y, i) = sin(x) / x;
      REAL(ky, i) = sin(x) / x; 
      x = (i * dx) + range[0];
      REAL(y, i + (N/2)) = sin(x) / x;
      REAL(ky, i + (N/2)) = sin(x) / x;
    }
  REAL(y,N/2) = sin(range[0])/range[0];
  REAL(y,N/2) = sin(range[0])/range[0];
  wavetable = gsl_fft_complex_wavetable_alloc (N);
  workspace = gsl_fft_complex_workspace_alloc (N);
  gsl_fft_complex_forward (ky, 1, N, wavetable, workspace);  // Taking fourier transform

  
  for (int i = 0; i < N/2; i++)      // Shifting the fourier plane
    {
      REAL(kys,i) = REAL(ky,i + (N/2));
      IMAG(kys,i) = IMAG(ky,i + (N/2));
      REAL(kys,i + (N/2)) = REAL(ky,i);
      IMAG(kys,i + (N/2)) = IMAG(ky,i);
    }
  
  FILE *fp;      // Saving to csv file
  int st=remove("P_3.csv");
  fp = fopen("P_3.csv", "w+");
  c_fac = dx/sqrt(2*M_PI);
  for(int i=0;i<N;i++)
    {
      k = 2*M_PI*((i - (N/2)) * dk);
      fprintf(fp, "%f,%f \n", k, c_fac*REAL(kys,i));
    }
  
  free(y); free(ky); free(kys); fclose(fp);
  gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);

  return(0);
}
