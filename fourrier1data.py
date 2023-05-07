import MyIntegration as mi
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

def FourierCoeff(expr, N, L, method, ftype = -1, d = 1, n = 100, pw = False,n_max = 1000):
    a0 = 0
    an = []
    bn = []

    if method == "trapezoidal" or method == "Trapezoidal":
        if ftype == 0:                      # even function
            a0 = mi.MyTrap_tol(expr,L,0,n_max,d)[0]/(L)

            for i in range(1,N+1):
                bn.append(0)
                func = "(" + expr + ")" + "*" + "cos({:}*(x*pi/{:}))".format(str(i), str(L))
                a_coeff = mi.MyTrap_tol(func,L,0,n_max,d)[0]
                an.append(2*a_coeff/L)
    
            if (pw == True):
                an = np.array(an)/2       # /2 because we require int-0 to L and not int-(-L) to L.
        
            an = np.array(an)
            bn = np.array(bn)

        elif ftype == 1:                       # odd function

            for i in range(1, N+1):
                an.append(0)
                func = "(" + expr + ")" + "*" + "sin({:}*(x*pi/{:}))".format(str(i), str(L))
                b_coeff = mi.MyTrap_tol(func,L,0,n_max,d)[0]
                bn.append(2*b_coeff/L)
            
            if (pw == True):
                a0 = mi.MyTrap_tol(expr,L,0,n_max,d)[0]/(2*L) # a0/2
                bn = np.array(bn)/2       # /2 because we require int-0 to L and not int-(-L) to L.

            an = np.array(an)
            bn = np.array(bn)
        
        elif ftype == -1:    # Neither Odd nor Even Function
            a0 = mi.MyTrap_tol(expr,L,-L,n_max,d)[0]/(2*L)    # a0/2

            if (pw == True):
                a0 = a0/2

                for i in range(1, N+1): 
                    func_cos = "(" + expr + ")" + "*" + "cos({:}*(x*pi/(2*{:})))".format(str(i), str(L))
                    a_coeff = mi.MyTrap_tol(func_cos,L,-L,n_max,d)[0]
                    an.append(a_coeff/(2*L))

                    func_sin = "(" + expr + ")" + "*" + "sin({:}*(x*pi/(2*{:})))".format(str(i), str(L))
                    b_coeff = mi.MyTrap_tol(func_sin,L,-L,n_max,d)[0]
                    bn.append(b_coeff/(2*L))
            
            else:
                for i in range(1, N+1): 
                    func_cos = "(" + expr + ")" + "*" + "cos({:}*(x*pi/{:}))".format(str(i), str(L))
                    a_coeff = mi.MyTrap_tol(func_cos,L,-L,n_max,d)[0]
                    an.append(a_coeff/L)

                    func_sin = "(" + expr + ")" + "*" + "sin({:}*(x*pi/{:}))".format(str(i), str(L))
                    b_coeff = mi.MyTrap_tol(func_sin,L,-L,n_max,d)[0]
                    bn.append(b_coeff/L)
            
            an = np.array(an)
            bn = np.array(bn)

    elif method == "simpson" or method == "Simpson":
        if ftype == 0:                      # even function
            a0 = mi.MySimp_tol(expr,L,0,n_max,d)[0]/(L)

            for i in range(1,N+1):
                bn.append(0)
                func = "(" + expr + ")" + "*" + "cos({:}*(x*pi/{:}))".format(str(i), str(L))
                a_coeff = mi.MySimp_tol(func,L,0,n_max,d)[0]
                an.append(2*a_coeff/L)
    
            if (pw == True):
                an = np.array(an)/2       # /2 because we require int-0 to L and not int-(-L) to L.
        
            an = np.array(an)
            bn = np.array(bn)

        elif ftype == 1:                       # odd function

            for i in range(1, N+1):
                an.append(0)
                func = "(" + expr + ")" + "*" + "sin({:}*(x*pi/{:}))".format(str(i), str(L))
                b_coeff = mi.MySimp_tol(func,L,0,n_max,d)[0]
                bn.append(2*b_coeff/L)
            
            if (pw == True):
                a0 = mi.MySimp_tol(expr,L,0,n_max,d)[0]/(2*L) # a0/2
                bn = np.array(bn)/2       # /2 because we require int-0 to L and not int-(-L) to L.

            an = np.array(an)
            bn = np.array(bn)
        
        elif ftype == -1:    # Neither Odd nor Even Function
            a0 = mi.MySimp_tol(expr,L,-L,n_max,d)[0]/(2*L)    # a0/2

            if (pw == True):
                a0 = a0/2

                for i in range(1, N+1): 
                    func_cos = "(" + expr + ")" + "*" + "cos({:}*(x*pi/(2*{:})))".format(str(i), str(L))
                    a_coeff = mi.MySimp_tol(func_cos,L,-L,n_max,d)[0]
                    an.append(a_coeff/(2*L))

                    func_sin = "(" + expr + ")" + "*" + "sin({:}*(x*pi/(2*{:})))".format(str(i), str(L))
                    b_coeff = mi.MySimp_tol(func_sin,L,-L,n_max,d)[0]
                    bn.append(b_coeff/(2*L))
            
            else:
                for i in range(1, N+1): 
                    func_cos = "(" + expr + ")" + "*" + "cos({:}*(x*pi/{:}))".format(str(i), str(L))
                    a_coeff = mi.MySimp_tol(func_cos,L,-L,n_max,d)[0]
                    an.append(a_coeff/L)

                    func_sin = "(" + expr + ")" + "*" + "sin({:}*(x*pi/{:}))".format(str(i), str(L))
                    b_coeff = mi.MySimp_tol(func_sin,L,-L,n_max,d)[0]
                    bn.append(b_coeff/L)
            
            an = np.array(an)
            bn = np.array(bn)

    elif method == "gauss" or method == "Gauss":
        if ftype == 0:                      # even function
            a0 = mi.MyLegQuadrature_tol(0,L,expr,2,d,n_max)[0]/(L)

            for i in range(1,N+1):
                bn.append(0)
                func = "(" + expr + ")" + "*" + "cos({:}*(x*pi/{:}))".format(str(i), str(L))
                a_coeff = mi.MyLegQuadrature_tol(0,L,expr,2,d,n_max)[0]
                an.append(2*a_coeff/L)
    
            if (pw == True):
                an = np.array(an)/2       # /2 because we require int-0 to L and not int-(-L) to L.
        
            an = np.array(an)
            bn = np.array(bn)

        elif ftype == 1:                       # odd function

            for i in range(1, N+1):
                an.append(0)
                func = "(" + expr + ")" + "*" + "sin({:}*(x*pi/{:}))".format(str(i), str(L))
                b_coeff = mi.MyLegQuadrature_tol(0,L,expr,2,d,n_max)[0]
                bn.append(2*b_coeff/L)
            
            if (pw == True):
                a0 = mi.MyLegQuadrature_tol(0,L,expr,2,d,n_max)[0]/(2*L) # a0/2
                bn = np.array(bn)/2       # /2 because we require int-0 to L and not int-(-L) to L.

            an = np.array(an)
            bn = np.array(bn)
        
        elif ftype == -1:    # Neither Odd nor Even Function
            a0 = mi.MyLegQuadrature_tol(-L,L,expr,2,d,n_max)[0]/(2*L)    # a0/2

            if (pw == True):
                a0 = a0/2

                for i in range(1, N+1): 
                    func_cos = "(" + expr + ")" + "*" + "cos({:}*(x*pi/(2*{:})))".format(str(i), str(L))
                    a_coeff = mi.MyLegQuadrature_tol(-L,L,expr,2,d,n_max)[0]
                    an.append(a_coeff/(2*L))

                    func_sin = "(" + expr + ")" + "*" + "sin({:}*(x*pi/(2*{:})))".format(str(i), str(L))
                    b_coeff = mi.MyLegQuadrature_tol(-L,L,expr,2,d,n_max)[0]
                    bn.append(b_coeff/(2*L))
            
            else:
                for i in range(1, N+1): 
                    func_cos = "(" + expr + ")" + "*" + "cos({:}*(x*pi/{:}))".format(str(i), str(L))
                    a_coeff = mi.MyLegQuadrature_tol(-L,L,expr,2,d,n_max)[0]
                    an.append(a_coeff/L)

                    func_sin = "(" + expr + ")" + "*" + "sin({:}*(x*pi/{:}))".format(str(i), str(L))
                    b_coeff = mi.MyLegQuadrature_tol(-L,L,expr,2,d,n_max)[0]
                    bn.append(b_coeff/L)
            
            an = np.array(an)
            bn = np.array(bn)

    return a0,an,bn

if __name__ == "__main__":
    print("\nName: ")

    n = [1, 2, 5, 10, 20]
    xx = np.linspace(-2.5, 2.5, 50)
    x = var("x")

    # # Function 1

    f1 = np.piecewise(xx, [ ((-3<xx) & (xx<-2)), ((-2<xx) & (xx<-1)), ((-1<xx) & (xx<0)), ((0<xx) & (xx<1)), ((1<xx) & (xx<2)), ((2<xx) & (xx<3))], [0, 1, 0, 1, 0, 1])
    dataf1 = []
    outf1 = []
    errf1 = []
    partialsum_f1 = []
    for i in n:
        a0, an, bn = FourierCoeff(expr = "1", N = i, L = 1,d = 4, method ='trapezoidal', ftype = 1, pw = True)
        partialsumf1 = a0
        for j in range(len(an)):
            partialsumf1 += bn[j]*sin(int(j+1)*np.pi*x) + an[j]*cos(int(j+1)*np.pi*x)
        partialsum_f1.append(lambdify(x, partialsumf1))
    
    plt.figure("Piecewise Function 1")
    plt.plot(xx, f1, label = "Function 1")
    for (p, q) in zip(partialsum_f1, n):
        plt.plot(xx, p(xx), marker = '.', linewidth = 0.5, label = "N = %d"%q)
        dataf1.append(p(xx))
        outf1.append(p(-0.5))
        outf1.append(p(0))
        outf1.append(p(0.5))
        errf1.append(np.abs(p(-0.5) - np.zeros(shape = [1, 50])))
        errf1.append(np.abs(p(0) - 0.5*np.ones(shape = [1, 50])))
        errf1.append(np.abs(p(0.5) - np.ones(shape = [1, 50])))


    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title("\nFourier Approximation")
    plt.legend()
    data_f1 = np.column_stack((xx, *dataf1))
    np.savetxt(r'F:name.txt', data_f1, header = "x, 1, 2, 5, 10, 20")
    data_out_f1 = np.column_stack((*outf1, *errf1))
    np.savetxt(r'F:name.txt', data_out_f1, header = "n=1 fr(-0.5), fr(0), fr(0.5), n=2 fr(-0.5), fr(0), fr(0.5), n=5 fr(-0.5), fr(0), fr(0.5), n=10 fr(-0.5), fr(0), fr(0.5), n=20 fr(-0.5), fr(0), fr(0.5), n=1 err(-0.5), err(0), err(0.5), n=2 err(-0.5), err(0), err(0.5), n=5 err(-0.5), err(0), err(0.5), n=10 err(-0.5), err(0), err(0.5), n=20 err(-0.5), err(0), err(0.5)")

    ## Fuction 2

    f2 = np.piecewise(xx, [ ((-3<xx) & (xx<-2.5)), ((-2.5<xx) & (xx<-1.5)), ((-1.5<xx) & (xx<-1)), ((-1<xx) & (xx<-0.5)), ((-0.5<xx) & (xx<0.5)), ((0.5<xx) & (xx<1)), ((1<xx) & (xx<1.5)), ((1.5<xx) & (xx<2.5)), ((2.5<xx) & (xx<3))], [0, 1, 0, 0, 1, 0, 0, 1, 0])
    dataf2 = []
    outf2 = []
    errf2 = []
    partialsum_f2 = []
    for i in n:
        a0, an, bn = FourierCoeff(expr = "1", N = i, L = 0.5,d = 4, method = 'trapezoidal', ftype = -1, pw = True)
        partialsumf2 = a0
        for j in range(len(an)):
            partialsumf2 += bn[j]*sin(int(j+1)*np.pi*x) + an[j]*cos(int(j+1)*np.pi*x)
        partialsum_f2.append(lambdify(x, partialsumf2))

    plt.figure("Piecewise Function 2")
    plt.plot(xx, f2, label = "Function 2")
    for (p, q) in zip(partialsum_f2, n):
        plt.plot(xx, p(xx), marker = '.', linewidth = 0.5, label = "N = %d"%q)
        dataf2.append(p(xx))
        outf2.append(p(-0.5))
        outf2.append(p(0))
        outf2.append(p(0.5))
        errf2.append(np.abs(p(-0.5) - np.zeros(shape = [1, 50])))
        errf2.append(np.abs(p(0) - 0.5*np.ones(shape = [1, 50])))
        errf2.append(np.abs(p(0.5) - np.ones(shape = [1, 50])))
    
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title("\nFourier Approximation of Piecewise Function (2)")
    plt.legend(loc = "lower right")
    data_f2 = np.column_stack((xx, *dataf2))
    np.savetxt(r'name.txt', data_f2, header = "x, 1, 2, 5, 10, 20")
    data_out_f2 = np.column_stack((*outf2, *errf2))
    np.savetxt(r'name.txt', data_out_f2, header = "n=1 fr(-0.5), fr(0), fr(0.5), n=2 fr(-0.5), fr(0), fr(0.5), n=5 fr(-0.5), fr(0), fr(0.5), n=10 fr(-0.5), fr(0), fr(0.5), n=20 fr(-0.5), fr(0), fr(0.5), n=1 err(-0.5), err(0), err(0.5), n=2 err(-0.5), err(0), err(0.5), n=5 err(-0.5), err(0), err(0.5), n=10 err(-0.5), err(0), err(0.5), n=20 err(-0.5), err(0), err(0.5)")

    # # Function 3

    f3 = np.piecewise(xx, [ ((-3<xx) & (xx<-2)), ((-2<xx) & (xx<-1)), ((-1<xx) & (xx<0)), ((0<xx) & (xx<1)), ((1<xx) & (xx<2)), ((2<xx) & (xx<3))], [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5])
    dataf3 = []
    outf3 = []
    errf3 = []    
    partialsum_f3 = []
    for i in n:
        a0, an, bn = FourierCoeff(expr = "0.5", N = i, L = 1,d = 4, method = 'trapezoidal', ftype = 1, pw = False)
        partialsumf3 = a0
        for j in range(len(an)):
            partialsumf3 += bn[j]*sin(int(j+1)*np.pi*x) + an[j]*cos(int(j+1)*np.pi*x)
        partialsum_f3.append(lambdify(x, partialsumf3))

    plt.figure("Piecewise Function 3")
    plt.plot(xx, f3, label = "Function 3")
    for (p, q) in zip(partialsum_f3, n):
        plt.plot(xx, p(xx), marker = '.', linewidth = 0.5, label = "N = %d"%q)
        dataf3.append(p(xx))
        outf3.append(p(-0.5))
        outf3.append(p(0))
        outf3.append(p(0.5))
        errf3.append(np.abs(p(-0.5) - np.zeros(shape = [1, 50])))
        errf3.append(np.abs(p(0) - 0.5*np.ones(shape = [1, 50])))
        errf3.append(np.abs(p(0.5) - np.ones(shape = [1, 50])))
    
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title("\nFourier Approximation of Piecewise Function (3)")
    plt.legend()
    data_f3 = np.column_stack((xx, *dataf3))
    np.savetxt(r'F:name.txt', data_f3, header = "x, 1, 2, 5, 10, 20")
    data_out_f3 = np.column_stack((*outf3, *errf3))
    np.savetxt(r'name.txt', data_out_f3, header = "n=1 fr(-0.5), fr(0), fr(0.5), n=2 fr(-0.5), fr(0), fr(0.5), n=5 fr(-0.5), fr(0), fr(0.5), n=10 fr(-0.5), fr(0), fr(0.5), n=20 fr(-0.5), fr(0), fr(0.5), n=1 err(-0.5), err(0), err(0.5), n=2 err(-0.5), err(0), err(0.5), n=5 err(-0.5), err(0), err(0.5), n=10 err(-0.5), err(0), err(0.5), n=20 err(-0.5), err(0), err(0.5)")

    plt.show()

