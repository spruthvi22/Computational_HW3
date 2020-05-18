 /*
* Author:Pruthvi Suryadevara
* Email: pruthvi.suryadevara@tifr.res.in
* Description: Finding Fourier Transform of Sinc function using FFTW
* Compile using gcc P_2.c -lfftw3 -lm -o P_2.out
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
  int N = 1024;
  double x, k, c_fac; 
  double range[] = {-100, 100};
  double dx = (range[1] - range[0]) / N;
  double dk = 1 / (N * dx);

  y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  ky = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  kys = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  y[0] = 1;

  for (int i = 1; i < N/2; i++)       // Defining shifted y 
    {
      x = (i * dx);
      y[i] = sin(x)/x;
      x = (i * dx) + range[0];
      y[i + (N/2)] = sin(x)/x;
    }
  y[N/2] = sin(range[0])/range[0];
  
  p = fftw_plan_dft_1d(N, y, ky, FFTW_FORWARD, FFTW_ESTIMATE);  // Taking fourier transform
  fftw_execute(p);
  
  c_fac = dx/sqrt(2*M_PI);
  
  for (int i = 0; i < N/2; i++)   // Shifting the fourier plane
    {
      kys[i] = ky[i + (N/2)];
      kys[i + (N/2)] = ky[i];
    }
  
  FILE *fp;
  int st=remove("P_2.csv");
  fp = fopen("P_2.csv", "w+"); 
  
  for(int i=0;i<N;i++)
    {
      k = 2*M_PI*((i - (N/2)) * dk);
      fprintf(fp, "%f,%f \n", k, c_fac*creal(kys[i]), creal(y[i]), cimag(y[i]));
    }
  
  fftw_destroy_plan(p);
  fftw_free(y); fftw_free(ky); fftw_free(kys); fclose(fp);

  return(0);
}
