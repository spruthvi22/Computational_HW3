 /*
* Author:Pruthvi Suryadevara
* Email: pruthvi.suryadevara@tifr.res.in
* Description: Finding Fourier transform of gaussian function
* Compile using gcc P_4.c -lfftw3 -lm -o P_4.out
* Note: Instead of multipiying the solution by phase factor we are defining y that is already shifted
* which mwke no difference as phase in fourier plane is shift in x plane
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

int main()
{
  fftw_complex *y, *ky, *kys;
  fftw_plan p;
  int N = 256;
  double x, c_factor; 
  double range[] = {-20, 20};
  double dx = (range[1] - range[0]) / (N);
  double dk = 1 / (N * dx);
  y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  ky = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  kys = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  double *act_sol = (double*) malloc(sizeof(double) * N);
  double *k = (double*) malloc(sizeof(double) * N);
     
  for (int i=0; i<N/2; i++)           // Defining shifted y gaussian function
    {
      x = (i*dx);
      y[i] = exp( -1 * x * x);
      x = (i * dx) + range[0];
      y[i +(N/2)] = exp( -1 * x * x);
    }
    
  p = fftw_plan_dft_1d(N, y, ky, FFTW_FORWARD, FFTW_ESTIMATE);  // Taking fourier transform
  fftw_execute(p); 

  for (int i = 0; i < N/2; i++)       // Shifting the fourier plane
    {
      kys[i] = ky[i + (N/2)];
      kys[i + (N/2)] = ky[i];
    }

  for (int i = 0; i < N; i++)        // Finding frequency and actual soluton
    {
      k[i] = 2*M_PI*((i - (N/2)) * dk);
      act_sol[i] = exp(-0.25 * k[i] * k[i]) / sqrt(2);
    }

  c_factor = dx/sqrt(2*M_PI);       // Coeffient for getting fourier transform form FFT
  FILE *fp;
  int st=remove("P_4.csv");
  fp = fopen("P_4.csv", "w+"); 

  for(int i=0;i<N;i++)              // Saving to file, Plotting done using Python
    {
      fprintf(fp, "%f,%f,%f \n", k[i], c_factor*creal(kys[i]), act_sol[i]);
    }
  
  fftw_destroy_plan(p);
  fftw_free(y); fftw_free(ky); free(k); fclose(fp);

  return(0);
}
