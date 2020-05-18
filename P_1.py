"""
Author: Pruthvi Suryadevara
Email:  pruthvi.suryadevara@tifr.res.in
Finding fourier transform of Sinc function
Note: Instead of multipiying the solution by phase factor we are defining y that is already shifted
which mwke no difference as phase in fourier plane is shift in x plane
 
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as ft


def f(x):                    # Defining Sinc Function
    y = np.zeros(len(x))
    for i in range(len(x)):
        if(x[i] != 0):
            y[i] = (np.sin(x[i]))/x[i]
        else:
            y[i] = 1
    return(y)


n = 1024                            # Taking number of points as power of 2
x = np.linspace(-100, 100, n+1)     # Defining x in the range
y = f(x[0:n])
y = ft.fftshift(y)                  # Defiing a symmetric matrix with y[0] is y(x=0)
ky = ft.fftshift(ft.fft(y))         # Fiding Fourier transform
k = 2*np.pi*ft.fftshift(ft.fftfreq(n, x[1]-x[0]))  # getting Frequency points   
c_factor = (x[1]-x[0])/np.sqrt(2*np.pi)    # constant factor to get Fourier transform form FFT
k_ac = np.linspace(k[0],k[len(k)-1],100)
a_sol=np.zeros(len(k_ac))
for i in range(len(k_ac)):
    if k_ac[i]>-1 and k_ac[i]<1:
        a_sol[i] = np.sqrt(np.pi/2)

plt.plot(k, c_factor*np.real(ky), 'b', label="Using numpy.fft")

p2 = np.genfromtxt('P_2.csv', delimiter=',')      # Loading solution by FFTW
p3 = np.genfromtxt('P_3.csv', delimiter=',')      # Loading solution form GSL
plt.plot(p2[:, 0], p2[:, 1], 'k', label="Using FFTW")
plt.plot(p3[:, 0], p3[:, 1], 'y', label="Using GSL")
plt.plot(k_ac, a_sol, '.r', label="Analytical solution")
plt.legend()
plt.xlabel("Frequency")
plt.ylabel("Fourier transform")
plt.title("Fourier tranform of Sinc Function")
plt.show()
