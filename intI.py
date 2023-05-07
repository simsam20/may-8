import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy as sp


def MyTrap(f,ul,ll,n):
    integrand = eval("lambda x :"+ f)
    nodal = [ll] ; func = []   # ARRAYS FOR CALCULATIONS
    # STEP SIZE
    h = (ul - ll)/n
    # NODAL POINTS
    for k in range(n):
        ll = ll + h
        nodal.append(ll)
    # FUNCTION VALUE AT NODAL POINT
    for j in nodal:
        func.append(integrand(j))
    # TRAPEZOIDAL INTEGRATION
    z2,inte = 0,0
    z1 = (func[0] + func[len(func)-1])/2
    for i in range(1,len(func)-1):
        z2 = z2 + func[i]
        inte = h * (z1 + z2)
    return inte
    
def MySimp(f,ul,ll,n):
    integrand = eval("lambda x :" + f)
    nodal = [ll] ; func = []   # ARRAYS FOR CALCULATIONS
    # STEP SIZE
    h = (ul - ll)/n
    # NODAL POINTS
    for k in range(n):
        ll = ll + h
        nodal.append(ll)
    # FUNCTION VALUE AT NODAL POINT
    for j in nodal:
        func.append(integrand(j))
    # SIMPSON INTEGRATION
    z2,z3,inte = 0,0,0
    if n%2 == 0:
        z1 = func[0] + func[len(func)-1]
        for i in range(1,len(func)-1):
            if i%2 == 0:
                z2 = z2 + func[i]
            else:
                z3 = z3 + func[i]
        inte = h/3 * (z1 + (4 * (z3)) + (2 * (z2)))
    else:
        print("Simpson Rule is not applicable since interval is not even.")
        inte = None
    return inte
    
def MyTrap_tol(f,ul,ll,n_max,d):
    tolerance = 0.1 * 10**(-d)
    count = [] ; I = []
    for i in range(2,n_max+4,2):
        I.append(MyTrap(f,ul,ll,i))
        count.append(i)
        if len(I) > 1:
            if abs((I[-1] - I[-2])/I[-1]) <= tolerance and i <= n_max:
                return [I[-1],count[-1]]
            elif i > n_max:
                print("The required tolerance could not be achieved with the maximum number of intervals.")
 
def MySimp_tol(f,ul,ll,n_max,d):
    tolerance = 10**(-d)
    count = [] ; I = []
    for i in range(2,n_max+4,2):
        I.append(MySimp(f,ul,ll,i))
        count.append(i)
        if len(I) > 1:
            if abs((I[-1] - I[-2])/I[-1]) <= tolerance and i <= n_max:
                return [I[-1],count[-1]]
            elif i > n_max:
                print("The required tolerance could not be achieved with the maximum number of intervals.")
                
def gauss_quad(a,b,N,f):
	x = var('x')
	expr = sympify(f)
	integrand = lambdify(x,expr)
	m = (b-a)/2
	c = (b+a)/2
	z = np.polynomial.legendre.leggauss(N)#z array is the array of points and weights
	return m*np.dot(z[1],integrand(m*z[0]+c))#Integral value

def MyLegQuadrature(a,b,N,f,m):
	h = (b-a)/m
	Int = []
	result = 0 
	for j in range(m+1):
		Int.append(a+j*h)
	for k in range(len(Int)-1):
		result += gauss_quad(Int[k],Int[k+1],N,f)
	return result

def MyLegQuadrature_tol(a,b,f,N,d,m_max):
    tolerance = 0.5*10**(-d)
    i = 1
    m = 1
    err = abs((MyLegQuadrature(a,b,N,f,1)-MyLegQuadrature(a,b,N,f,2))/MyLegQuadrature(a,b,N,f,2))
    while err>=tolerance and m<=m_max:
        i += 1
        memo1= MyLegQuadrature(a,b,N,f,m)
        m = 2**i
        memo2 = MyLegQuadrature(a,b,N,f,m)
        err = abs((memo2-memo1)/memo2)
        if m>m_max:
            print("the sub intervals exceeded the given maximum limit")
            return None 
    result = [MyLegQuadrature(a,b,N,f,2**i) , m ]
    return result

def MyLaguQuad(expression, n):
    func = eval("lambda x :"+ expression)

    xvals,weights  = np.polynomial.laguerre.laggauss(n)
    w_fx = [weights[i]*func(xvals[i]) for i in range(len(weights))]
    result = sum(w_fx)

    return result

def MyHermiteQuad(expression, n):
    func = eval("lambda x :"+ expression)

    xvals,weights  = np.polynomial.hermite.hermgauss(n)
    w_fx = [weights[i]*func(xvals[i]) for i in range(len(weights))]
    result = sum(w_fx)

    return result
    

                


