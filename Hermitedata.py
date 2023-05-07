import numpy as np
import matplotlib.pyplot as plt
import MyIntegration as mi
import pandas as pd



def integral_simp(I1_exact,I2_exact):
    i = 2
    I1_simp = []
    I2_simp = []
    count1_simp = []
    count2_simp = []
    I1_exact_S = []
    I2_exact_S = []

    while True:
        if abs(I1_exact - mi.MySimp("exp(-x**2)/(1+x**2)",1000,-1000,i)) <= 10**(-2):
            I1_simp.append(mi.MySimp("exp(-x**2)/(1+x**2)",1000,-1000,i))
            I1_exact_S.append(I1_exact)
            count1_simp.append(i)
            break
        else:
            I1_simp.append(mi.MySimp("exp(-x**2)/(1+x**2)",1000,-1000,i))
            I1_exact_S.append(I1_exact)
            count1_simp.append(i)
            i += 2

    i = 2

    while True:
        if abs(I2_exact - mi.MySimp("1/(1+x**2)",1000,-1000,i)) <= 10**(-2):
            I2_simp.append(mi.MySimp("1/(1+x**2)",1000,-1000,i))
            I2_exact_S.append(I2_exact)
            count2_simp.append(i)
            break
        else:
            I2_simp.append(mi.MySimp("1/(1+x**2)",1000,-1000,i))
            I2_exact_S.append(I2_exact)
            count2_simp.append(i)
            i += 2

    return I1_simp,I2_simp,count1_simp,count2_simp,I1_exact_S,I2_exact_S

def graph(I1,I2,n,I1_simp,I2_simp,count1_simp,count2_simp,I1_exact_S,I2_exact_S):
    fig1,ax1 = plt.subplots(1, 2)
    fig2,ax2 = plt.subplots(1, 2)
    ax1[0].plot(n,I1,label = "MyHermiteQuad")
    ax1[0].plot(n,I1_exact_LL,label = "Analytic Value")
    ax1[1].plot(count1_simp,I1_simp,label = "MySimp")
    ax1[1].plot(count1_simp,I1_exact_S,label = "Analytic Value")
    ax2[0].plot(n,I2,label = "MyHermiteQuad")
    ax2[0].plot(n,I2_exact_LL,label = "Analytic Value")
    ax2[1].plot(count2_simp,I2_simp,label = "MySimp")
    ax2[1].plot(count2_simp,I2_exact_S,label = "Analytic Value")
    for i in range(2):
        if i == 0:
            ax1[i].set(xlabel = "Nodal Points (n)",ylabel = "Value of Integration (I)",title = "Gauss Hermite Quadrature")
            ax2[i].set(xlabel = "Nodal Points (n)",ylabel = "Value of Integration (I)",title = "Gauss Hermite Quadrature")
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
    Exact=[1.7724538509,0,0.886226925,0,1.329340388,0,3.32335097044,0,11.63172839,0]
    
    for count in range(len(Exact)):
        f = input("\nEnter Function: ")
        func.append(f)
        if count < (len(Exact) - 1):
            ans = input("Do you want to enter more function (Y/N) ?\t")
            if ans == "N" or ans == "n":
                break

    for j,m in zip(func,Exact):
        for k in range(2,6,2):
            print("\nValue of integration of",j,"for n =",k,"is: ",mi.MyHermiteQuad(j,k))
        print("\nExact Value of integration of",j,"is: ",m)
        print("---------------------------------------------------------------")

    
    # PART B II

    I1 = []
    I2 = []
    I1_exact = 1.343293421646735
    I2_exact = 3.141592653589793
    I1_exact_LL = []
    I2_exact_LL = []
    n = []
    
    for i in range(2,130,2):
        n.append(i)
        I1_exact_LL.append(I1_exact)
        I2_exact_LL.append(I2_exact)
        i1 = mi.MyHermiteQuad("1/(1+x**2)",i)
        I1.append(i1)
        i2 = mi.MyHermiteQuad("exp(x**2)/(1+x**2)",i)
        I2.append(i2)

    data1 =  np.column_stack([n,I1,I2])
    file1 = np.savetxt("quad-herm-1221.txt",data1,header = ("n,I1,I2"))

    df1 =  pd.DataFrame({"n": n, "I1": I1, "I2": I2})
    print("\nGAUSS HERMITE QUADRATURE:\n",df1)

    # PART B III & IV

    I1_simp,I2_simp,count1_simp,count2_simp,I1_exact_S,I2_exact_S = integral_simp(I1_exact,I2_exact)

    df2 =  pd.DataFrame({"n": count1_simp, "I1": I1_simp})
    print("\nTOLERNACE LIMIT = 10**(-2)")
    print("\nSIMPSON FOR INTEGRAL 1:\n",df2)
    data2 =  np.column_stack([count1_simp,I1_simp])
    file2 = np.savetxt("Simpson-Integral_1(H)-1221.txt",data2,header = ("n,I1"))

    print("\n-----------------------------------------------------------------")

    df3 =  pd.DataFrame({"n": count2_simp, "I2": I2_simp})
    print("\nTOLERNACE LIMIT = 10**(-2)")
    print("\nSIMPSON FOR INTEGRAL 1:\n",df3)
    data3 =  np.column_stack([count2_simp,I2_simp])
    file3 = np.savetxt("Simpson-Integral_2(H)-1221.txt",data3,header = ("n,I2"))

    graph(I1,I2,n,I1_simp,I2_simp,count1_simp,count2_simp,I1_exact_S,I2_exact_S)

