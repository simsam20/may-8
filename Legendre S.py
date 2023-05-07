import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre

def MyLegendre(n, x):
    y = np.zeros(len(x))
    for i in range(n+1):
        y += legendre(i)(x) * coeff[i]
    return y

def find_coeff(f, n):
    coeff = np.zeros(n+1)
    for i in range(n+1):
        f_legendre = lambda x: f(x) * legendre(i)(x)
        coeff[i] = 2 / (1.0 * (i + 1)) * np.trapz(f_legendre(x), x)
    return coeff

def f1(x):
    return 2 + 3*x + 4*x**4

def f2(x):
    return np.cos(x) * np.sin(x)

n1 = 6
n2 = 8
x = np.linspace(-2, 2, 1000)

coeff = find_coeff(f1, n1)
y1 = MyLegendre(n1, x)

coeff = find_coeff(f2, n2)
y2 = MyLegendre(n2, x)

fig, axs = plt.subplots(2, 1, figsize=(8, 8))

axs[0].set_title('f(x) = 2 + 3x + 4x^4')
axs[0].plot(x, f1(x), label='f(x)')
for i in range(1, n1+1):
    axs[0].plot(x, coeff[i]*legendre(i)(x), label=f'n = {i}')
axs[0].legend()

axs[1].set_title('f(x) = cos(x) * sin(x)')
axs[1].plot(x, f2(x), label='f(x)')
for i in range(2, n2+1, 2):
    axs[1].plot(x, coeff[i]*legendre(i)(x), label=f'n = {i}')
axs[1].legend()

plt.show()
