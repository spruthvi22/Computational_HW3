"""
Author: Pruthvi Suryadevara
Email:  pruthvi.suryadevara@tifr.res.in
Time diffreence between DFT and FFT

"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as ft
import time

# Taking 8 points that are powere of 2
N_range = np.round(np.geomspace(2**5, 2**12, 8)).astype(int)
fft_time = []
dft_time = []
dft2_time = []
for N in N_range:
    np.random.seed(0)
    x = np.linspace(0, 10, N)
    Y = np.random.rand(N)          # Creating a random y of size n, to FFT
    strt_time = time.time()
    Yk_fft = ft.fft(Y)             # Perofomring FFT and recording time
    end_time = time.time()
    fft_time.append(end_time - strt_time)
    # Performing DFT and recording time, two times are recorded as
    # one which includes generation of M matrix and other without
    # as some libraries have M matrices stored in storage to improve time
    strt_time1 = time.time()
    M = np.zeros([N, N], dtype=np.complex)
    i = np.arange(N)                         
    k = np.reshape(i, (N, 1))
    M = np.exp(2j * i * k * np.pi / N)
    strt_time2 = time.time()
    Yk_dft = np.dot(M, Y)
    end_time1 = time.time()
    dft_time.append(end_time1 - strt_time1)
    dft2_time.append(end_time1 - strt_time2)


plt.subplot(2, 1, 1, title="time of dft vs fft")
plt.plot(N_range, np.array(fft_time), 'r', label="fft")
plt.plot(N_range, np.array(dft2_time), 'k', label="dft without M matrix")
plt.plot(N_range, np.array(dft_time), 'b', label="dft with M matrix")
plt.legend()
plt.xlabel("Number of points (N)")
plt.ylabel("Time of execution (sec)")

plt.subplot(2, 1, 2, title="Fourier tranform using dft and fft")
plt.plot(x, np.real(Yk_fft), 'r', label="fft")
plt.plot(x, np.real(Yk_dft), 'b', label="dft")
plt.legend()
plt.xlabel("k")
plt.ylabel("fourier transform")
plt.subplots_adjust(hspace=0.5)
plt.show()
