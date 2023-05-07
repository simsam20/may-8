#20068567002
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd

def MySinSeries(x,n):
    sin_cal = 0
    for i in range(n):
        sin_cal = sin_cal + (-1)**i * x**(2*i + 1)/math.factorial(2*i + 1)
    return sin_cal

def MyCosSeries(x,n):
    cos_cal = 0
    for i in range(n):
        cos_cal = cos_cal + (-1)**i * x**(2*i)/ math.factorial(2*i)
    return cos_cal

def tol(x):
    sigd = int(input("Number of Significant Digits: "))
    tol = 10**-(sigd)
    arr = []
    arrsin = []
    for i in xx:
        n = 1
        if (tol >= np.abs((MySinSeries(i, n+1) - MySinSeries(i, n))/(MySinSeries(i, n)))):        
            arr.append(n)
            arrsin.append(MySinSeries(i, n))
        else:
            while(tol < np.abs((MySinSeries(i, n+1) - MySinSeries(i, n))/(MySinSeries(i, n)))):   
                n += 1
            arr.append(n)
            arrsin.append(MySinSeries(i, n))
    return(arr, tol)
        
def graph(y_s,y_c,angle,inb_s,inb_c,y0_c,y0_s,n,sine_i,cosine_i,my_sin,inbuilt_sin,xx):
    fig,ax =  plt.subplots()
    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots(2,1)
    fig3,ax3 = plt.subplots()
    for i in range(len(y_s)):
        ax.plot(angle,y_s[i],label = n_a[i])
        ax1.plot(angle,y_c[i],label = n_a[i])
    ax.plot(angle,inb_s,label = "Inbuilt Sine Funtion")
    ax1.plot(angle,inb_c,label = "Inbuilt Cosine Function")
    ax2[0].plot(n,y0_s,label = "My Sine Function")
    ax2[0].scatter(n,y0_s,marker = ".")
    ax2[0].plot(n,sine_i,label = "Inbuilt Sine")
    ax2[1].plot(n,y0_c,label = "My Cosine Function")
    ax2[1].scatter(n,y0_c,marker = ".")
    ax2[1].plot(n,cosine_i,label = "Inbuilt Cosine")
    ax3.plot(xx,my_sin,label = "My Sine Function")
    ax3.scatter(xx,my_sin,marker = ".")
    ax3.plot(xx,inbuilt_sin,label = "Inbuilt Sine Function")
    ax.legend()
    ax1.legend()
    ax2[0].legend()
    ax2[1].legend()
    ax3.legend()
    ax.set(ylim = (-10,10),title = "Sine Series",xlabel = "Angle (x)",ylabel = "sin(x)")
    ax1.set(ylim = (-10,10),title = "Cosine Series",xlabel = "Angle (x)",ylabel = "cos(x)")
    ax2[0].set(xlabel = "Number of Terms",ylabel = "sin($x_o$)")
    ax2[1].set(xlabel = "Number of Terms",ylabel = "cos($x_o$)")
    ax3.set(xlabel = "Angle (x)",ylabel = "sin(x)",title = "Precision Graph")
    ax.grid(ls = "--")
    ax1.grid(ls = "--")
    ax2[0].grid(ls = "--")
    ax2[1].grid(ls = "--")
    ax3.grid(ls = "--")
    fig2.suptitle("$x_o$ = $\pi/4$")
    plt.show()

if __name__ == "__main__":
    n_a = [1,2,5,10,20] # PART A
    n_b = 20 ; n =[] ; sine_i = [] ;cosine_i = [] # PART B & C
    
    angle = np.arange(-2*np.pi,2*np.pi,0.1)
    y_s = [] ; y_c = [] ; y0_s = [] ; y0_c = []
    
    inb_s = np.sin(angle)
    inb_c = np.cos(angle)

    for j in range(len(n_a)): # SOL PART A
        out_s = MySinSeries(angle,n_a[j])   
        y_s.append(out_s)
        out_c = MyCosSeries(angle,n_a[j])
        y_c.append(out_c)
        
    for j in range(2,n_b+1): # SOL PART B & C
        n.append(j)
        sine_i.append(np.sin(np.pi/4)) ; cosine_i.append(np.cos(np.pi/4))
        sine = MySinSeries((np.pi/4),j)
        y0_s.append(sine)
        cosine = MyCosSeries((np.pi/4),j)
        y0_c.append(cosine)
        
    xx = np.arange(0,np.pi,np.pi/8) # Question 2
    index = 0
    tolarr, tol1 = tol(xx)
    inbuilt_sin = np.sin(xx)
    my_sin = []
    for p in xx:
        my_sin.append(MySinSeries(p, tolarr[index]))
        index += 1
    print("\n Tolerance = ",tol1)
    
    dict = {'x': xx,
            'sin(x) calc': my_sin,
            'Number of terms (n)': tolarr, 
            'sin(x) inbuilt': inbuilt_sin} 
    df = pd.DataFrame(dict)
    print(df)
        
    graph(y_s,y_c,angle,inb_s,inb_c,y0_c,y0_s,n,sine_i,cosine_i,my_sin,inbuilt_sin,xx)
        


