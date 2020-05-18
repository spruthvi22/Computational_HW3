"""
Author: Pruthvi Suryadevara
Email:  pruthvi.suryadevara@tifr.res.in
Finding convolution of box function with itself

"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as ft


def f(x):                             # Defining Step Function
    y = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] < 1 and x[i] > -1:
            y[i] = 1
        else:
            y[i] = 0
    return(y)


n = 512
x = np.linspace(-10, 10, n+1)[0:n]
y = f(x)
k = ft.fft(y, norm='ortho')                         # Taking FFT
k = k**2
conv = (x[1]-x[0])*np.sqrt(n)*ft.fftshift(ft.ifft(k, norm='ortho'))

plt.plot(x, np.real(conv), 'k', label="Convolution")
plt.plot(x, y, 'b', label="Step Function")
plt.xlabel("x")
plt.ylabel("Convolution")
plt.title("Convoluion of step function with itself")
plt.legend()
plt.show()
