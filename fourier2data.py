from test import FourierCoeff
import matplotlib.pyplot as plt
import numpy as np
from sympy.abc import x
from sympy import sin, cos, lambdify

if __name__ == '__main__':
    print("\nName:")

    n = [1, 2, 5, 10, 20]
    xx = np.linspace(-np.pi, np.pi, 50)
    
    # Odd Extension

    f_odd = xx
    outodd = []
    errodd = []
    partialsum_odd = []

    for i in n:
        a0, an, bn = FourierCoeff(expr = "x", N = i, L = np.pi,d = 4, method ='simpson', ftype = 1, pw = False)
        partialsumodd = a0
        for j in range(len(an)):
            partialsumodd += bn[j]*sin(int(j+1)*np.pi*x/np.pi) + an[j]*cos(int(j+1)*np.pi*x/np.pi)
        partialsum_odd.append(lambdify(x, partialsumodd))
    
    plt.figure("x odd")
    plt.plot(xx, f_odd, label = "Odd Extension: f(x) = x")
    for (p, q) in zip(partialsum_odd, n):
        plt.plot(xx, p(xx), marker = '.', linewidth = 0.5, label = "N = %d"%q)
        outodd.append(p(0))
        outodd.append(p(np.pi/2))
        outodd.append(p(np.pi))
        errodd.append(np.abs(p(0) - np.zeros(shape = [1, 50])))
        errodd.append(np.abs(p(np.pi/2) - np.pi/2*np.ones(shape = [1, 50])))
        errodd.append(np.abs(p(np.pi) - np.pi*np.ones(shape = [1, 50])))

    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title("")
    plt.legend()
    data_odd = np.column_stack((*outodd, *errodd))
    np.savetxt(r'.txt', data_odd, header = "n=1 fr(-0.5), fr(0), fr(0.5), n=2 fr(-0.5), fr(0), fr(0.5), n=5 fr(-0.5), fr(0), fr(0.5), n=10 fr(-0.5), fr(0), fr(0.5), n=20 fr(-0.5), fr(0), fr(0.5), n=1 err(-0.5), err(0), err(0.5), n=2 err(-0.5), err(0), err(0.5), n=5 err(-0.5), err(0), err(0.5), n=10 err(-0.5), err(0), err(0.5), n=20 err(-0.5), err(0), err(0.5)")

    # Even Extension

    f_even = np.abs(xx)
    outeven = []
    erreven = []
    partialsum_even = []

    for i in n:
        a0, an, bn = FourierCoeff(expr = "x", N = i, L = np.pi,d = 4, method ='trapezoidal', ftype = 0, pw = False)
        partialsumeven = a0
        for j in range(len(an)):
            partialsumeven += bn[j]*sin(int(j+1)*np.pi*x/np.pi) + an[j]*cos(int(j+1)*np.pi*x/np.pi)
        partialsum_even.append(lambdify(x, partialsumeven))
    
    plt.figure("x even")
    plt.plot(xx, f_even, label = "Even Extension: f(x) = x")
    for (p, q) in zip(partialsum_even, n):
        plt.plot(xx, p(xx), marker = '.', linewidth = 0.5, label = "N = %d"%q)
        outeven.append(p(0))
        outeven.append(p(np.pi/2))
        outeven.append(p(np.pi))
        erreven.append(np.abs(p(0) - np.zeros(shape = [1, 50])))
        erreven.append(np.abs(p(np.pi/2) - np.pi/2*np.ones(shape = [1, 50])))
        erreven.append(np.abs(p(np.pi) - np.pi*np.ones(shape = [1, 50])))

    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title("\nFourier Approximation of Even Extension of f(x) = x")
    plt.legend()
    data_even = np.column_stack((*outodd, *errodd))
    np.savetxt(r'', data_even, header = "n=1 fr(-0.5), fr(0), fr(0.5), n=2 fr(-0.5), fr(0), fr(0.5), n=5 fr(-0.5), fr(0), fr(0.5), n=10 fr(-0.5), fr(0), fr(0.5), n=20 fr(-0.5), fr(0), fr(0.5), n=1 err(-0.5), err(0), err(0.5), n=2 err(-0.5), err(0), err(0.5), n=5 err(-0.5), err(0), err(0.5), n=10 err(-0.5), err(0), err(0.5), n=20 err(-0.5), err(0), err(0.5)")

    plt.show()
