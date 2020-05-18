"""
Author: Pruthvi Suryadevara
Email:  pruthvi.suryadevara@tifr.res.in
Finding fourier transform of 2D Gaussian Function

"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as ft


def f(x):
    return(np.exp(-1 * x))


N = 256
x = np.linspace(-50, 50, N+1)[0:N]
y = x
X, Y = np.meshgrid(x, y)
z = f(X**2 + Y**2)
plt.show()
z = ft.fftshift(z)
zk = ft.fft2(z)
kx = 2*np.pi*ft.fftshift(ft.fftfreq(len(x), x[1]-x[0]))
ky = kx
kX, kY = np.meshgrid(kx, ky)
zk = ft.fftshift(zk)
c_factor = (x[1] - x[0])*(y[1] - y[0])/(2*np.pi)
act_sol = 0.5*np.exp(-0.25*(kX**2 + kY**2))
print(np.max(act_sol))
print(np.max(c_factor*np.real(zk)))
plt.subplot(1, 2, 1, title="Using DFT")
plt.imshow(c_factor*np.real(zk), extent=[kx[0], kx[len(kx)-1], ky[0], ky[len(ky)-1]], aspect='auto')
plt.colorbar()
plt.subplot(1, 2, 2, title="Actual Soluction")
plt.imshow(act_sol, extent=[kx[0], kx[len(kx)-1], ky[0], ky[len(ky)-1]], aspect='auto')
plt.colorbar()
plt.show() 
