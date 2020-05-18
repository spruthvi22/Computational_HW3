"""
Author: Pruthvi Suryadevara
Email:  pruthvi.suryadevara@tifr.res.in
Plotting Fourier tranform of gaussian function

"""

import numpy as np
import matplotlib.pyplot as plt


p4 = np.genfromtxt('P_4.csv', delimiter=',')    # Importing Data from file and printing
plt.plot(p4[:, 0], p4[:, 1], 'k', label="Using FFTW")
plt.plot(p4[:, 0], p4[:, 2], '.r', label="Analutical solution")
plt.legend()
plt.xlabel("Frequency")
plt.ylabel("Fourier transform")
plt.title("Fourier Transform of Gaussian fucntion")
plt.show()
