"""
Author: Pruthvi Suryadevara
Email:  pruthvi.suryadevara@tifr.res.in
Finding Power Spectrum Density of a noise

"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as ft
import scipy.stats as st


y = np.loadtxt("noise.txt")
ky = ft.fftshift(ft.fft(y))                         # Taking FFT
k = 2*np.pi*ft.fftshift(ft.fftfreq(len(y), 1))      # Finding Frequency
plt.subplot(2, 2, 1, title="Noise")
plt.plot(y)
plt.subplot(2, 2, 2, title="DFT")
plt.plot(k, np.real(ky), label="Real part")
plt.plot(k, np.imag(ky), label="Imaginary part")
plt.legend()
plt.subplot(2, 2, 3, title="Periodogram")
Pky = (np.abs(ky)**2)/(len(y)*2*np.pi)              # Finding Periodogram and plotting
plt.plot(k, Pky)
plt.subplot(2, 2, 4, title="Binned Periodogram")    # Binning Periodogram
ky_bin, k_be, binnumber = st.binned_statistic(k, Pky, bins=10)
k_bins = (k_be[0:len(k_be)-1]+k_be[1:len(k_be)])/2
plt.bar(k_bins, ky_bin, width=k_be[1]-k_be[0])
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.show()
