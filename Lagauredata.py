''' import numpy as np


My Function for integrating a function using trapezoidal method that works for both – fixed number of intervals and fixed tolerance.
Input Parameters -  function (f) , a - lower limit , b - upper limit , n0 - Number of panels , key1(Bool)= False(for finding Integral for only 1 value of n(intervals)) True(for finding integral for more than 1 value of n) , N_max=maximum value of n , key2(Bool)= True (for tolerance) , tol=tolerance

def MyTrap(f,a,b,n0,key1=True,N_max=None,key2=False,tol=None):
    w=0  
    h=(b-a)/n0         #step size
    S=0.5*(f(a)+f(b))          
    for i in range(1,n0):
        S+= f(a+i*h)
    Integral = S * h
    if key1== True:
        pass
    else:
        return Integral      #returning the integral if key is set to false
    n_a=[n0]       #creating array for values of n(intervals)
    I=[Integral]       #creating array to store values of integral
    n0=2*n0          #doubling subintervals
    while n0<=N_max:
        h=(b-a)/n0
        r=0
        for i in range(1,n0):        #calculating the value of function at certain points to avoid repeated calculations 
              if i%2 != 0:
                 r+=f(a+i*h)
        r=r*h
        I.append((I[-1]/2)+(r))     
        n_a.append(n0)
        if key2==True:
            if I[-1]<=0.1e-25:
                err=abs(I[-1]-I[-2])
            else:
               err=abs((I[-1]-I[-2])/I[-1])  #calculation of relative error
            if err<=tol:
                w=1         
                break
            else:
                pass
        else:
            pass
        n0=2*n0        
    if key2==True:    #printing the message if key2 is true    
        if w==0:
            s=("N_max reached without achieving required tolerance")
        elif w==1:
             s="Given tolerance achieved with",n_a[-1],"intervals"
        return I[-1],n_a[-1],s       #returning integral,number of intervals and message
    else:
        return I[-1],n_a[-1]          #returning integral,number of intervals
        

My Function for integrating a function using simpson 1/3 method that works for both – fixed number of intervals and fixed tolerance.
Input Parameters -  function (f) , a - lower limit , b - upper limit , n0 - Number of panels , key1(Bool)= False(for finding Integral for only 1 value of n(intervals)) True(for finding integral for more than 1 value of n) , N_max=maximum value of n , key2(Bool)= True (for tolerance) , tol=tolerance


def MySimp(f,a,b,n0,key1=True,N_max=None,key2=False,tol=None):
    if n0%2 ==0:
        pass
    else :
        return "Number of intervals must be even"     #this works for even number of sub intervals
    w=0
    S_a=[];T_a=[];I_a=[]
    h=(b-a)/n0      #step size
    S = f(a) + f(b)
    T=0
    for i in range(1,n0):       
        if i%2 == 0:
            S = S + 2 * f(a + i*h)
        else:
            T = T + (2 * f(a + i*h))/3
    S=S/3
    Integral =h*(S+2*T)
    if key1== True:
        pass
    else:
        return Integral   #returning the integral if key is set to false
    S_a.append(S);T_a.append(T);I_a.append(Integral);n_a=[n0]        #creating array for values of n(intervals),Integral
    n0=2*n0  #doubling subintervals
    while n0<=N_max:
        h=(b-a)/n0
        T=0
        for i in range(1,n0):      #calculating the value of function at certain points to avoid repeated calculations 
            if i%2 != 0:
                T+=(2*f(a+i*h))/3
        
        S=S_a[-1]+T_a[-1]
        Integral =h*(S+2*T)
        S_a.append(S);T_a.append(T);I_a.append(Integral);n_a.append(n0)
        if key2==True:
            if I_a[-1]<=0.1e-25:
                err=abs(I_a[-1]-I_a[-2])
            else:
                err=abs((I_a[-1]-I_a[-2])/I_a[-1])           #calculation of relative error
            if err<=tol:
                w=1
                break
            else:
                pass
        else:
            pass
        n0=2*n0
    if key2==True:      #printing the message if key2 is true  
        if w==0:
            s=("N_max reached without achieving required tolerance")
        elif w==1:
             s="Given tolerance achieved with",n_a[-1],"intervals"
        return I_a[-1],n_a[-1],s     #returning integral,number of intervals and message
    else:
        return I_a[-1],n_a[-1]      #returning integral,number of intervals
    
My Function for integrating a function using Legendre Gauss method that works for any order and any sub intervals and tolerance is optional parameter
Input Parameters -  function (f) , a - lower limit , b - upper limit ,  n - order of gauss legendre , m - Number of sub intervals , key(Bool)= False (for finding Integral for only a certain value of m and n ) True(for finding integral for more than 1 value of n,m upto certain tolerance) , m_max=maximum value of m  , tol=tolerance



def MyLegQuadrature(f,a,b,n,m,key=False,tol=None,m_max=None):
  def gs1(f,a,b,n,m0):           #subfunction which returns the value of integral for specific n and m
      x,w = np.polynomial.legendre.leggauss(n)
      h=(b-a)/m
      s=0
      for i in range(0,m):
          r=0
          for x1,w1 in zip(x,w):
               r+=w1*f((((a+(i+1)*h)-(a+i*h))/2)*x1+((a+i*h)+(a+(i+1)*h))/2)
          r= (((a+(i+1)*h)-(a+i*h)) /2 )*r
          s+=r
      return s
  Integral=gs1(f,a,b,n,m)
  if key==True:
     pass
  else:
     
     return Integral     #returning thr value of integral if key is false
  w=0
  m_a=[m]
  I=[Integral]
  m=2*m
  while m<=m_max:
    I.append(gs1(f,a,b,n,m))
    m_a.append(m)
    if I[-1]<=0.1e-25:
        err=abs(I[-1]-I[-2])
    else:
        err=abs((I[-1]-I[-2])/I[-1])
    if err<=tol:
        w=1
        break
    else:
        pass
    
    m=2*m
  if w==0:
      s=("m_max reached without achieving required tolerance")
  elif w==1:
      s="Given tolerance achieved with",m_a[-1],"sub-intervals"
  return [I[-1],m_a[-1],s]        #returns integral,number of subintervals and message

def MyLaguQuad(f, n):
    Integral=0
    xi,wi=np.polynomial.laguerre.laggauss(n)
    for (Xi,Wi) in zip (xi,wi):
        Integral+=Wi*f(Xi)
    return Integral

def MyHermiteQuad(f, n):
    Integral=0
    xi,wi=np.polynomial.hermite.hermgauss(n)
    for (Xi,Wi) in zip (xi,wi):
        Integral+=Wi*f(Xi)
    return Integral'''

import numpy as np
import matplotlib.pyplot as plt
import MyIntegration as mi
import pandas as pd

# Name: Ishmeet Singh, 2020PHY1221
# Partner Name: Sarthak Jain, 2020PHY1201"

def integral_simp(I1_exact,I2_exact):
    i = 2
    I1_simp = []
    I2_simp = []
    count1_simp = []
    count2_simp = []
    I1_exact_S = []
    I2_exact_S = []

    while True:
        if abs(I1_exact - mi.MySimp("exp(-x)/(1+x**2)",1000,0,i)) <= 10**(-2):
            I1_simp.append(mi.MySimp("exp(-x)/(1+x**2)",1000,0,i))
            I1_exact_S.append(I1_exact)
            count1_simp.append(i)
            break
        else:
            I1_simp.append(mi.MySimp("exp(-x)/(1+x**2)",1000,0,i))
            I1_exact_S.append(I1_exact)
            count1_simp.append(i)
            i += 2

    i = 2

    while True:
        if abs(I2_exact - mi.MySimp("1/(1+x**2)",1000,0,i)) <= 10**(-2):
            I2_simp.append(mi.MySimp("1/(1+x**2)",1000,0,i))
            I2_exact_S.append(I2_exact)
            count2_simp.append(i)
            break
        else:
            I2_simp.append(mi.MySimp("1/(1+x**2)",1000,0,i))
            I2_exact_S.append(I2_exact)
            count2_simp.append(i)
            i += 2

    return I1_simp,I2_simp,count1_simp,count2_simp,I1_exact_S,I2_exact_S

def graph(I1,I2,n,I1_simp,I2_simp,count1_simp,count2_simp,I1_exact_S,I2_exact_S):
    fig1,ax1 = plt.subplots(1, 2)
    fig2,ax2 = plt.subplots(1, 2)
    ax1[0].plot(n,I1,label = "MyLagQuad")
    ax1[0].plot(n,I1_exact_LL,label = "Analytic Value")
    ax1[1].plot(count1_simp,I1_simp,label = "MySimp")
    ax1[1].plot(count1_simp,I1_exact_S,label = "Analytic Value")
    ax2[0].plot(n,I2,label = "MyLagQuad")
    ax2[0].plot(n,I2_exact_LL,label = "Analytic Value")
    ax2[1].plot(count2_simp,I2_simp,label = "MySimp")
    ax2[1].plot(count2_simp,I2_exact_S,label = "Analytic Value")
    for i in range(2):
        if i == 0:
            ax1[i].set(xlabel = "Nodal Points (n)",ylabel = "Value of Integration (I)",title = "Gauss Laguerre Quadrature")
            ax2[i].set(xlabel = "Nodal Points (n)",ylabel = "Value of Integration (I)",title = "Gauss Laguerre Quadrature")
        elif i == 1:
            ax1[i].set(xlabel = "Nodal Points (n)",ylabel = "Value of Integration (I)",title = "Simpson 1/3 Method")
            ax2[i].set(xlabel = "Nodal Points (n)",ylabel = "Value of Integration (I)",title = "Simpson 1/3 Method")
        ax1[i].grid(ls = "--")
        ax2[i].grid(ls = "--")
        ax1[i].legend()
        ax2[i].legend()
    fig1.suptitle("INTEGRAL 1")
    fig2.suptitle("INTEGRAL 2")
    plt.show()
        
if __name__ == "__main__":
    print("\nName: Ishmeet Singh\tRoll No. : 2020PHY1221\nPartner: Sarthak Jain\tRoll No.: 2020PHY1201\n")

    # PART B I

    count = 0
    func = []
    Exact=[1,1,2,6,24,120,720,5040,40320,362880]
    
    for count in range(len(Exact)):
        f = input("\nEnter Function: ")
        func.append(f)
        if count < (len(Exact) - 1):
            ans = input("Do you want to enter more function (Y/N) ?\t")
            if ans == "N" or ans == "n":
                break

    for j,m in zip(func,Exact):
        for k in range(2,6,2):
            print("\nValue of integration of",j,"for n =",k,"is: ",mi.MyLaguQuad(j,k))
        print("\nExact Value of integration of",j,"is: ",m)
        print("---------------------------------------------------------------")

    
    # PART B II

    I1 = []
    I2 = []
    I1_exact = 0.621449624235813
    I2_exact = 1.570796326794897
    I1_exact_LL = []
    I2_exact_LL = []
    n = []
    
    for i in range(2,130,2):
        n.append(i)
        I1_exact_LL.append(I1_exact)
        I2_exact_LL.append(I2_exact)
        i1 = mi.MyLaguQuad("1/(1+x**2)",i)
        I1.append(i1)
        i2 = mi.MyLaguQuad("exp(x)/(1+x**2)",i)
        I2.append(i2)

    data1 =  np.column_stack([n,I1,I2])
    file1 = np.savetxt("quad-lag-1221.txt",data1,header = ("n,I1,I2"))

    df1 =  pd.DataFrame({"n": n, "I1": I1, "I2": I2})
    print("\nGAUSS LAGUERRE QUADRATURE:\n",df1)

    # PART B III & IV

    I1_simp,I2_simp,count1_simp,count2_simp,I1_exact_S,I2_exact_S = integral_simp(I1_exact,I2_exact)

    df2 =  pd.DataFrame({"n": count1_simp, "I1": I1_simp})
    print("\nTOLERNACE LIMIT = 10**(-2)")
    print("\nSIMPSON FOR INTEGRAL 1:\n",df2)
    data2 =  np.column_stack([count1_simp,I1_simp])
    file2 = np.savetxt("Simpson-Integral_1-1221.txt",data2,header = ("n,I1"))

    print("\n-----------------------------------------------------------------")

    df3 =  pd.DataFrame({"n": count2_simp, "I1": I2_simp})
    print("\nTOLERNACE LIMIT = 10**(-2)")
    print("\nSIMPSON FOR INTEGRAL 1:\n",df3)
    data3 =  np.column_stack([count2_simp,I2_simp])
    file3 = np.savetxt("Simpson-Integral_2-1221.txt",data3,header = ("n,I2"))

    graph(I1,I2,n,I1_simp,I2_simp,count1_simp,count2_simp,I1_exact_S,I2_exact_S)