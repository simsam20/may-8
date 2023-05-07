from random import gauss
from sympy.abc import x
import numpy as np
from numpy.polynomial.legendre import leggauss
from sympy import simplify, poly, degree
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.hermite import hermgauss

def getval(expr, **kwargs):
    func = eval("lambda x:"+ expr)      # Lambdifying the function

    for key in kwargs.keys():

        if ((key == "x") or (key == "xarr") or (key == "xint")):
            xarr = list(kwargs.get(key))
            a = xarr[0]
            b = xarr[-1]
            n = len(xarr)

            if (("tol" in kwargs.keys()) and (kwargs.get("tol") != None)):
                tol = kwargs.get("tol")
            else:
                tol = 0.5*10**(-8)    # Default tol

            break
        
        elif ((key == "a") or (key == "b")):
            a = kwargs.get("a")
            b = kwargs.get("b")
            n = 100                   # Default n 
            tol = 0.5*10e-8           # Default tol      

            if (("n" in kwargs.keys()) and (kwargs.get("n") != None)):
                n = kwargs.get("n")

            if (("tol" in kwargs.keys()) and (kwargs.get("tol") != None)):
                tol = kwargs.get("tol")  

            break
        
        else:
            raise(KeyError("Bounds not found! Provide Lower and Upper Bounds or the Complete Interval to proceed"))
            
    return(func, a, b, n, tol)

def trapztol(tol, y1, func, h1, a, n):
    n_max = 2**20
    trap1 = h1*sum(y1)

    n = 2*n
    h2 = h1/2
    y2 = np.zeros(n)
    
    for p in range(1, n, 2):
        y2[p] = func(a+p*h2)
    for q in range(int (n/2)):
        y2[2*q] = y1[q]
    trap2 = h2*sum(y2)

    if (tol >= np.abs((trap2 - trap1)/trap2)):
        result_n = n/2

    else:
        while (tol < np.abs((trap2 - trap1)/trap2)):
            if(n <= n_max):
                y3 = np.zeros(2*n)
                h2 = h2/2
                trap1 = trap2

                for p in range(1, 2*n, 2):
                    y3[p] = func(a+p*h2)
                for q in range(n):
                    y3[2*q] = y2[q]
                trap2 = h2*sum(y3)
                y1 = y2
                y2 = y3
                n *= 2 

                result_n = n/2

            else:
                print("Tolerance could not be achieved with n(max) = ", n_max)
                result_n = n/2
                break

    return (result_n, y1)

def simpstol(tol, y_1, func, h1, a, n):
    y1 = y_1
    n_max = 2**20
    h2 = h1/2
    y2 = np.zeros(2*n)
    result1 = h1*sum(y1)
    result_n = n

    for i in range (0, n, 2):
        y2[2*i] = y1[i]
    for j in range (1, n, 2):
        y2[2*j] = y1[j]/2
    for k in range(1, 2*n, 2):
        y2[k] = 4*func(a + k*h2)/3
    result2 = h2*sum(y2)

    n = 2*n
    
    if (tol >= np.abs((result2 - result1)/result2)):
        result_n = n/2

    else:
        while (tol < np.abs((result2 - result1)/result2)):
            if (n <= n_max):
                y3 = np.zeros(2*n)
                h2 = h2/2
                result1 = result2
                for i in range (0, n, 2):
                    y3[2*i] = y2[i]
                for j in range (1, n, 2):
                    y3[2*j] = y2[j]/2
                for k in range(1, 2*n, 2):
                    y3[k] = 4*func(a + k*h2)/3
                result2 = h2*sum(y3)
                y1 = y2
                y2 = y3
                n *= 2

                result_n = n/2
            else:
                print("Tolerance could not be achieved with n(max) = ", n_max)
                result_n = n/2
                break

    return (result_n, y1)

def gausstol (func, tol, a, b, m, rootx, weightx, result1):
    m_max = 2**15
    m = 2*m
    h2 = (b-a)/m
    a2 = a
    result2 = 0
    result_m = m
    for i in range(m):
        b2 = a2 + h2
        new_xvals = (h2*rootx + (b2+a2))/2
        w_fx = [weightx[j]*func(new_xvals[j]) for j in range(len(weightx))]
        result2 += h2*sum(w_fx)/2
        a2 = b2

    if (tol >= np.abs((result2 - result1)/result1)):
        result_m = m/2

    else:
        while(tol < np.abs((result2 - result1)/result1)):
            if(m <= m_max):
                a3 = a
                result3 = 0
                h2 = h2/2
                for i in range(2*m):
                    b3 = a3 + h2
                    new_xvals = (h2*rootx + (b3+a3))/2
                    w_fx = [weightx[j]*func(new_xvals[j]) for j in range(len(weightx))]
                    result3 += h2*sum(w_fx)/2
                    a3 = b3
                m *= 2
                result1 = result2
                result2 = result3
                result_m = m/2
    
            else:
                print("Tolerance could not be achieved with m(max) = ", m_max)
                result_m = m/2
                break
    
    return(result_m, result1)

def trapz(expr, **kwargs):
    func, a, b, n, tol = getval(expr, **kwargs)
    h1 = (b-a)/n
    y1 = [(func(a)+func(b))/2]

    for i in range(1, n):
        y1.append(func(a+i*h1))         # Value of f(x) at the nodal points
    tol_n, y_final = trapztol(tol, y1, func, h1, a, n)
    trap = ((b-a)/tol_n) * sum(y_final)

    return(trap, tol_n)

def simps(expr,**kwargs):
    func, a, b, n, tol = getval(expr, **kwargs)
    if (n%2 != 0):
        raise(ValueError("n must be even!"))
        exit

    h = (b-a)/(n)
    y = [(func(a)+func(b))/3]
    
    for i in range(1, n): 
        if(i%2 == 0):
            y.append(2*func(a+i*h)/3)  # y at Even Nodes (at even indexes)         
        elif(i%2 == 1):
            y.append(4*func(a+i*h)/3)  # y at Odd Nodes (at odd indexes)
    tol_n, y_final = simpstol(tol, y, func, h, a, n)
    simp = ((b-a)/tol_n)*sum(y_final)
    
    return(simp, tol_n)

def gaussquad(expression, n = None, m = None, **kwargs):
    func, ll, ul, *other_params = getval(expression, **kwargs)
    a = ll  # Lower Limit
    b = ul  # Upper Limit

    if (n == None):
        if (simplify(expression).is_polynomial() == False):
            raise(KeyError("Expression is not a Polynomial. To Approximate the integral, enter n."))
        else:
            expr = poly(expression, x, domain = 'ZZ')
            xvals, weights = leggauss(degree(expr, gen = x))
            xvals = ((b-a)*xvals + (b+a))/2
            w_fx = [weights[i]*func(xvals[i]) for i in range(len(weights))]
        result = (b-a)*sum(w_fx)/2
        tol_m = None

    elif (n >= 1):
        xvals, weights = leggauss(n)
        tol = other_params[1]
        if (m != None):
            h = float((b - a)/m)
            result = 0
            for i in range(m):
                b = a + h
                new_xvals = (h*xvals + (b+a))/2
                w_fx = [weights[j]*func(new_xvals[j]) for j in range(len(weights))]
                result += h*sum(w_fx)/2
                a = b
            tol_m, result_final = gausstol(func, tol, ll, ul, m, xvals, weights, result)
            result = round(result_final, 14)        # np.sum often messes with floating point numbers

        else:
            xvals = ((b-a)*xvals + (b+a))/2
            w_fx = [weights[i]*func(xvals[i]) for i in range(len(weights))]
            result = (b-a)*sum(w_fx)/2
            tol_m = None

    return(result, tol_m)

def gausslag(expression, n):
    func = eval("lambda x:" + expression)

    xvals, weights = laggauss(n)
    w_fx = [weights[i]*func(xvals[i]) for i in range(len(weights))]
    result = sum(w_fx)

    return(result)

def gausshermite(expression, n):
    func = eval("lambda x:" + expression)

    xvals, weights = hermgauss(n)
    w_fx = [weights[i]*func(xvals[i]) for i in range(len(weights))]
    result = sum(w_fx)

    return(result)

if __name__ == '__main__':
    repeat = 'y'
    while(repeat == 'y'):
        function = str(input("\nEnter the Function f(x): "))
        # xx = np.linspace(1, 2, 100)
        # t, nn = gaussquad(expression = function, a = -1, b = 1, n = 5, m = 1, tol = 1)
        # t, nn = trapz(function, a = 0, b = 1, n = 4, tol = 1)
        # t, nn = simps(function, a = -100000, b = 100000, n = 2000000, tol = 1)
        t = gausshermite(function, n = 10)
        print("Integral is = ", t)
        # print("Tolerance = 10e-5 reached with n = ", nn)
        repeat = str(input("\nAnother Integration? Enter y/n : "))

