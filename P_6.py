"""
Author: Pruthvi Suryadevara
Email:  pruthvi.suryadevara@tifr.res.in
Finding fourier transform of Constant function

"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as ft


def FFT(x):
    #A recursive implementation FFT
    x = np.asarray(x, dtype=float)
    N = x.shape[0]
    if N % 2 > 0:
        raise ValueError("size of x must be a power of 2")
    elif N <= 2:  # Cutoff at N=2, returning known FFT
        return np.array([(x[0]+x[1]), (x[0]-x[1])])
    else:
        X_even = FFT(x[0:N:2])                    # Finding FFT of even indices
        X_odd = FFT(x[1:N:2])                     # Finding FFT of off indices
        factor = np.exp(-2j * np.pi * np.arange(N) / N)    
        return np.concatenate([X_even + (factor[0:int(N / 2)] * X_odd), X_even+ (factor[int(N / 2):N] * X_odd)])

x = np.linspace(-10, 10, 512)
y = np.ones(len(x))
ky = ft.fftshift(ft.fft(y))          # FFT using library function
k = 2*np.pi*ft.fftshift(ft.fftfreq(len(x), x[1]-x[0]))
c_factor = (x[1]-x[0])/np.sqrt(2*np.pi)   # Constant factor for findinf fourier tranform form FFT
ky2 = FFT(y)                         # FFT using defined function
ky2 = ft.fftshift(ky2)

plt.plot(k, c_factor*np.real(ky), 'b', label="Using library function")
plt.plot(k, c_factor*np.real(ky2), '.r', label="Using Defined function")
plt.legend()
plt.xlabel("Frequency")
plt.ylabel("Fourier transform")
plt.title("Fourier tranform Constant function")
plt.show()
