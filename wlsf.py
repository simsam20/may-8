import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def Mylsf(x,y):
    n = len(x)
    # SLOPE
    slope = (n*np.sum(x*y) - (np.sum(x)*np.sum(y)))/(n*np.sum((x**2)) - (np.sum(x))**2)
    
    #INTERCEPT
    intercept = (((np.sum(np.power(x,2)))*np.sum(y)) - (np.sum(x*y)*np.sum(x)))/((n*np.sum(np.power(x,2))) - (np.power(np.sum(x),2)))
    
    # Y CALCULATED
    y_cal = []
    for i in range(n):
        y_cal.append(slope*x[i] + intercept)
    y_cal = np.array(y_cal)
    
    # SUM OF RESIDUAL
    res = []
    for j in range(n):
        res.append(y[j] - y_cal[j])
    res = np.array(res)
    res_s = np.sum(res)
    
    # SQUARE OF SUM oF RESIDUAL
    res_1 = []
    for j in range(n):
        res_1.append((y[j] - y_cal[j])**2)
    res_1 = np.array(res_1)
    res_ss = np.sum(res_1)
    
    # STANDARD ERROR IN SLOPE
    s = np.sqrt(res_ss/(n-2))
    SS_xx = ((np.sum(np.power(x,2))) - (np.power(np.sum(x),2)/n))
    Serr_slope = s/(np.sqrt(SS_xx))
    
    # STANDARD ERROR IN INTERCEPT
    Serr_intercept = (Serr_slope * ((np.sqrt(np.sum(np.power(x,2))))/n))
    
    # COFFICIENT OF DETERMINATION
    diff_l = []
    y_mean = np.mean(y)
    for k in range(n):
        diff_l.append((y[k] - y_mean)**2)
    diff_l = np.array(diff_l)
    tss = np.sum(diff_l)
    R_sqr = (tss - res_ss)/tss # cofficient of determination
    
    # CORRELATION COFFICIENT
    cc = np.sqrt(R_sqr) # cofficient of correlation
    
    return slope,intercept,y_cal,res_s,res_ss,Serr_slope,Serr_intercept,cc

def MyWlsf(x,y,w):
    n = len (x)
    wxy = []
    wx = []
    wy = []
    wx2 = []
    for i in range (n):
       wxy.append(w[i]* x[i]* y[i])
       wx.append(w[i]* x[i])
       wy.append(w[i]* y[i])
       wx2.append(w[i]* x[i]* x[i])
    
    # SLOPE
    slope_w = ( np.sum (w)* np.sum ( wxy ) - np.sum ( wx )* np.sum ( wy )) /( np.sum (w) *np.sum ( wx2 ) - np.sum( wx )* np.sum ( wx ))
    slope_w = np.array(slope_w)

    # INTERCEPT
    intercept_w = (np.sum ( wx2 )* np.sum ( wy ) - np.sum( wx )* np.sum ( wxy )) /( np.sum( w)* np.sum ( wx2 ) - np.sum ( wx )* np.sum ( wx ))
    intercept_w = np.array(intercept_w)

    # Y CALCULATED
    y_cal_w = []
    for i in range(n):
        y_cal_w.append(slope_w*x[i] + intercept_w)
    y_cal_w = np.array(y_cal_w)

    # SUM OF RESIDUAL
    res_w = []
    for j in range(n):
        res_w.append(y[j] - y_cal_w[j])
    res_w = np.array(res_w)
    res_s_w = np.sum(res_w)
    
    # SQUARE OF SUM oF RESIDUAL
    res_1_w = []
    for j in range(n):
        res_1_w.append((y[j] - y_cal_w[j])**2)
    res_1_w = np.array(res_1_w)
    res_ss_w = np.sum(res_1_w)

    # ERROR IN SLOPE
    slope_err_wls = np.sqrt(( np.sum (w)) /( np.sum(w )* np.sum( wx2 ) - np.sum( wx )* np.sum( wx )))
    
    # ERROR IN INTERCEPT
    intercept_err_wls = np.sqrt(( np.sum ( wx2 )) /( np.sum(w)* np.sum( wx2 ) - np.sum( wx )* np.sum( wx ) ))

    # COFFICIENT OF DETERMINATION
    diff_l = []
    y_mean = np.mean(y)
    for k in range(n):
        diff_l.append((y[k] - y_mean)**2)
    diff_l = np.array(diff_l)
    tss = np.sum(diff_l)
    R_sqr_wls = (tss - res_ss_w)/tss # cofficient of determination

    # CORRELATION COFFICIENT
    cc_wls = np.sqrt(R_sqr_wls) # cofficient of correlation

    return slope_w,intercept_w,y_cal_w,res_s_w,res_ss_w,slope_err_wls,intercept_err_wls,R_sqr_wls,cc_wls


if __name__ == "__main__":

    print()

    x_vals = pd.read_csv(r"name of csv", usecols = [1])
    x_vals = (x_vals.to_numpy()).flatten()
    df = pd.read_csv(r"name of csv", usecols = range(2, 12))
    df = df.to_numpy()
    

    y_mean = np.array([])
    y_std_error = np.array([])
    for i in range(len(df[0])):
        mean = np.mean(df[i])
        y_mean = np.append(y_mean, mean)
        var = np.var(df[i])
        std_error = (4*mean**2)*var/len(df[0])
        y_std_error = np.append(y_std_error, std_error)

    x = x_vals.reshape(-1, 1)
    y = (y_mean**2).reshape(-1, 1)
    w = 1/y_std_error
   
    slope,intercept,y_cal,res_s,res_ss,Serr_slope,Serr_intercept,cc = Mylsf(x,y)
    slope_w,intercept_w,y_cal_w,res_s_w,res_ss_w,slope_err_wls,intercept_err_wls,R_sqr_wls,cc_wls = MyWlsf(x,y,w)

    data = np.column_stack([x,y,w])
    np.savetxt (" 1221. txt ",data, header = "xi , yi , wi")

    print("\nFITTING PARAMETERS (OLSF):")
    print("\nSlope: ",slope)
    print("\nError in Slope: ",Serr_slope)
    print("\nIntercept: ",intercept)
    print("\nError in Intercept: ",Serr_intercept)
    print("\nSum of residulas: ",res_s)
    print("\nSum of square of residulas: ",res_ss)
    print("\nCofficient of Correlation: ",cc)

    print("\n-----------------------------------------------------------------\n")

    print("\nFITTING PARAMETERS (WLSF):")
    print("\nSlope: ",slope_w)
    print("\nError in Slope: ",slope_err_wls)
    print("\nIntercept: ",intercept_w)
    print("\nError in Intercept: ",intercept_err_wls)
    print("\nSum of residulas: ",res_s_w)
    print("\nSum of square of residulas: ",res_ss_w)
    print("\nCofficient of Correlation: ",cc_wls)

    k = (4*np.pi**2)/slope_w
    m = intercept_w*k/(4*np.pi**2)
    error_k = slope_err_wls*k/slope_w
    error_m = (intercept_err_wls/intercept_w + error_k/k)*m

    print("\n------------------------------------------------------------------\n")

    print("\nValue of k: ",k)
    print("\nValue of m: ",m)
    print("\nError in k: ",error_k)
    print("\nError in m: ",error_m)

    # Plots and Scatters
    plt.scatter(x, y, marker = 'o')    
    plt.plot(x, y_cal, linestyle = 'dashed', linewidth = 1, label = "OLS Fitted Line",c = "red")
    plt.plot(x, y_cal_w, linewidth = 1, label = "WLS Fitted Line",c = "green")
    plt.title("")       #NAME OF EXPERIMENT
    plt.ylabel("")         #Y LABEL
    plt.xlabel("")                           #X LABEL
    plt.grid(ls = "--")

    plt.legend()
    plt.show()


    
    

