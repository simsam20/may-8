
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve

def finite_difference(a,b,alpha,beta,px,qx,rx,n,bc):
    A = np.zeros((n+1,n+1),float)
    B = np.zeros(n+1,float)

    if  bc == "R":
        h = (b-a)/n

        p = [px(a + i*h) for i in range(n+1)]
        q = [qx(a + i*h) for i in range(n+1)]
        r = [rx(a + i*h) for i in range(n+1)]
        
        u = [-1 + h/2*(p[i]) for i in range(1,n)]
        d = [2 + h**2*q[i] for i in range(1,n)]
        l = [-1 - h/2*(p[i]) for i in range(1,n)]
        b_i = [-h**2*r[i] for i in range(1,n)]

        a_11 =  (2 + (h**2)*q[0]) + 2*h*(-1 - h/2*(p[0]))*(alpha[0]/alpha[1])
        a_12 = -2

        a_n1n1 = (2 + (h**2)*q[-1]) - 2*h*(-1 + h/2*(p[-1]))*(beta[0]/beta[1])
        a_n1n = -2

        b_1 = (-h**2)*r[0] + 2*h*(-1 - h/2*(p[0]))*(alpha[2]/alpha[1])
        b_n1 = (-h**2)*r[-1] - 2*h*(-1 + h/2*(p[-1]))*(beta[2]/beta[1])

        A[0,0] = a_11
        A[0,1] = a_12
        A[-1,-1] = a_n1n1
        A[-1,-2] = a_n1n

        B[0] = b_1
        B[-1] = b_n1

        for i in range(1,n):
            B[i] = b_i[i-1]
            for j in range(n+1):
                if i == j:
                    A[i,j] = d[j-1]
                elif i + 1 == j:
                    A[i,j] = u[j-2]
                elif i == j + 1:
                    A[i,j] = l[j]

    elif bc == "N":
        h = (b-a)/n

        p = [px(a + i*h) for i in range(n+1)]
        q = [qx(a + i*h) for i in range(n+1)]
        r = [rx(a + i*h) for i in range(n+1)]
        
        u = [-1 + h/2*(p[i]) for i in range(1,n)]
        d = [2 + h**2*q[i] for i in range(1,n)]
        l = [-1 - h/2*(p[i]) for i in range(1,n)]
        b_i = [-h**2*r[i] for i in range(1,n)]

        a_11 =  (2 + (h**2)*q[0])
        a_12 = 0

        a_n1n1 = (2 + (h**2)*q[-1])
        a_n1n = 0

        b_1 = (-h**2)*r[0] + 2*h*(-1 - h/2*(p[0]))*(alpha[2]/alpha[1])
        b_n1 = (-h**2)*r[-1] - 2*h*(-1 + h/2*(p[-1]))*(beta[2]/beta[1])

        A[0,0] = a_11
        A[0,1] = a_12
        A[-1,-1] = a_n1n1
        A[-1,-2] = a_n1n

        B[0] = b_1
        B[-1] = b_n1

        for i in range(1,n):
            B[i] = b_i[i-1]
            for j in range(n+1):
                if i == j:
                    A[i,j] = d[j-1]
                elif i + 1 == j:
                    A[i,j] = u[j-2]
                elif i == j + 1:
                    A[i,j] = l[j]

    elif bc == "D":
        h = (b-a)/n

        p = [px(a + i*h) for i in range(n+1)]
        q = [qx(a + i*h) for i in range(n+1)]
        r = [rx(a + i*h) for i in range(n+1)]
        
        u = [-1 + h/2*(p[i]) for i in range(1,n)]
        d = [2 + h**2*q[i] for i in range(1,n)]
        l = [-1 - h/2*(p[i]) for i in range(1,n)]
        b_i = [-h**2*r[i] for i in range(1,n)]

        a_11 = 1
        a_12 = 0

        a_n1n1 = 1
        a_n1n = 0

        b_1 = (alpha[2])
        b_n1 = (beta[2])

        A[0,0] = a_11
        A[0,1] = a_12
        A[-1,-1] = a_n1n1
        A[-1,-2] = a_n1n

        B[0] = b_1
        B[-1] = b_n1

        for i in range(1,n):
            B[i] = b_i[i-1]
            for j in range(n+1):
                if i == j:
                    A[i,j] = d[j-1]
                elif i + 1 == j:
                    A[i,j] = u[j-2]
                elif i == j + 1:
                    A[i,j] = l[j]

    Y = solve(A,B)

    return A,B,Y


def exact(x_array,sol):
    if sol == 1:
        sol_e = (3/8)*(np.sin(x_array))-np.cos(x_array)-(1/8)*np.sin(3*x_array)
    elif sol == 2:
        sol_e = np.sin(np.pi*x_array)
    elif sol == 3:
        sol_e = (np.sin(np.pi*x_array))**2
    return sol_e                

if __name__ == "__main__":
    a = 0 ; b = 1
    n = 10
    alpha = [1,0,0]
    beta = [1,0,0]
    x_array = np.linspace(a,b,n+1)

    px  = lambda x : 0
    qx  = lambda x : np.pi**2
    rx  = lambda x : -2*np.pi**2 * np.sin(np.pi*x)

    A,B,Y = finite_difference(a,b,alpha,beta,px,qx,rx,n,"D")
    sol_e = exact(x_array,2)
    
    print(A)
    print(B)
    print("\n",Y)
    print("\n",sol_e)

    plt.plot(x_array,Y)
    plt.plot(x_array,sol_e,ls = "--")
    plt.show()

    a = 0 ; b = np.pi/2
    n = 10
    alpha = [1,1,-1]
    beta = [0,1,1]
    x_array = np.linspace(a,b,n+1)

    px  = lambda x : 0
    qx  = lambda x : -1
    rx  = lambda x : np.sin(3*x)

    A,B,Y = finite_difference(a,b,alpha,beta,px,qx,rx,n,"R")
    sol_e = exact(x_array,1)
    
    print(A)
    print(B)
    print("\n",Y)
    print("\n",sol_e)

    plt.plot(x_array,Y)
    plt.plot(x_array,sol_e,ls = "--")
    plt.show()

    a = 0 ; b = 1
    n = 10
    alpha = [0,1,-1]
    beta = [0,1,2*np.sin(1)]
    x_array = np.linspace(a,b,n+1)

    px  = lambda x : 0
    qx  = lambda x : -x
    rx  = lambda x : (3 - x - x**2 - x**3)*(np.sin(x)) + 4*x*(np.cos(x)) 

    A,B,Y = finite_difference(a,b,alpha,beta,px,qx,rx,n,"N")
    sol_e = exact(x_array,3)
    
    print(A)
    print(B)
    print("\n",Y)
    print("\n",sol_e)

    plt.plot(x_array,Y)
    plt.plot(x_array,sol_e,ls = "--")
    plt.show()
    
    
    

    
    