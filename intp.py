import numpy as np

'''
My Function for integrating a function using trapezoidal method that works for both – fixed number of intervals and fixed tolerance.
Input Parameters -  function (f) , a - lower limit , b - upper limit , n0 - Number of panels , key1(Bool)= False(for finding Integral for only 1 value of n(intervals)) True(for finding integral for more than 1 value of n) , N_max=maximum value of n , key2(Bool)= True (for tolerance) , tol=tolerance
'''
def MyTrap(f,args):
    a,b,n0,key1,N_max,key2,tol=args
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
            if I[-1]<=1e-20:
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
        
'''
My Function for integrating a function using simpson 1/3 method that works for both – fixed number of intervals and fixed tolerance.
Input Parameters -  function (f) , a - lower limit , b - upper limit , n0 - Number of panels , key1(Bool)= False(for finding Integral for only 1 value of n(intervals)) True(for finding integral for more than 1 value of n) , N_max=maximum value of n , key2(Bool)= True (for tolerance) , tol=tolerance
'''
#def MySimp(f,a,b,n0,key1=True,N_max=None,key2=False,tol=None):
def MySimp(f,args):
    a,b,n0,key1,N_max,key2,tol=args
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
               if I_a[-1]<=1e-20:
                   err=abs(I_a[-1]-I_a[-2])
               else:
                   err=abs((I_a[-1]-I_a[-2])/I_a[-1])  #calculation of relative error
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
    
'''
My Function for integrating a function using Legendre Gauss method that works for any order and any sub intervals and tolerance is optional parameter
Input Parameters -  function (f) , a - lower limit , b - upper limit ,  n - order of gauss legendre , m - Number of sub intervals , key(Bool)= False (for finding Integral for only a certain value of m and n ) True(for finding integral for more than 1 value of n,m upto certain tolerance) , m_max=maximum value of m  , tol=tolerance
'''


#def MyLegQuadrature(f,a,b,n,m,key=False,tol=None,m_max=None):
def MyLegQuadrature(f,args):
  a,b,n,m,key,tol,m_max=args
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
    if I[-1]<=1e-20:
        err=abs(I[-1]-I[-2])
    else:
        err=abs((I[-1]-I[-2])/I[-1])  #calculation of relative error    
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
  return [I[-1],m_a[-1],s]        #returns integral,number of 
