import numpy as np
import matplotlib.pyplot as plt
# y''' + 3y' = 0

def f(yin, xval):
    y1, y2 = yin
    dydt = y2
    dy2dt = np.sin(3*xval)-y1
    eqn = np.array([dydt, dy2dt])
    return eqn

def elr(func, y0, a, b, n):
    h = (b - a) / n
    yarr = []
    xarr = []
    x_axis = np.linspace(a, b, n)
    yin = y0 
    xin = [a]
    for i in range(n):
        yarr.append(yin)
        xarr.append(xin)
        yn = [e1 + h * ele for (e1, ele) in zip(yin, func(yin, x_axis[i]))]
        yin = yn
        xin = [e1 + h * (i+1) for e1 in xin]
    yarr = np.array(yarr).reshape(-1, len(yin))

    res_matrix = []
    for i in range(yarr.shape[1]):
        y_i = list(yarr[:, i:i+1].flatten())
        res_matrix.append(y_i)
    return(res_matrix)

def rk2(func, y0, a, b, n):
    h = (b - a) / n
    yarr = []
    yin = y0
    xin = [a]
    x_axis = np.linspace(a, b, n)
    for i in range(n):
        yarr.append(yin)
        k1 = [h * ele for ele in func(yin, x_axis[i])]
        yn = [e1 + e2/2 for (e1, e2) in zip(yin, k1)]
        xn = [e1 + h/2 for e1 in xin]
        k2 = [h * ele for ele in func(yn, x_axis[i])]
        yn = [e1 + e2/2 for (e1, e2) in zip(yin, k2)]
        yf = [iny + e1 for (iny, e1) in zip(yin, k2)]
        yin = yf
        xin = [e1 + h/2 for e1 in xn]
    yarr=np.array(yarr).reshape(-1, len(yin))

    res_matrix = []
    for i in range(yarr.shape[1]):
        y_i = list(yarr[:, i:i+1].flatten())
        res_matrix.append(y_i)
    return(res_matrix)

def rk4(func, y0, a, b, n):
    h = (b - a) / n
    yarr = []
    yin = y0
    x_axis = np.linspace(a, b, n+1)
    for i in range(n):
        yarr.append(yin)
        k1 = [h * ele for ele in func(yin, x_axis[i])]
        yn = [e1 + e2/2 for (e1, e2) in zip(yin, k1)]
        k2 = [h * ele for ele in func(yn, x_axis[i]+h/2)]
        yn = [e1 + e2/2 for (e1, e2) in zip(yin, k2)]
        k3 = [h * ele for ele in func(yn, x_axis[i]+h/2)]
        yn = [e1 + e2 for (e1, e2) in zip(yin, k3)]
        k4 = [h * ele for ele in func(yn, x_axis[i]+h)]
        yf = [ini_y + (e1 + 2 * (e2 + e3) + e4) / 6 for (ini_y,e1,e2,e3,e4) in zip(yin, k1, k2, k3, k4)]
        yin = yf
    yarr = np.array(yarr).reshape(-1, len(yin))

    res_matrix = []
    for i in range(yarr.shape[1]):
        y_i = list(yarr[:, i:i+1].flatten())
        res_matrix.append(y_i)
    return(res_matrix)

if __name__ == '__main__':
    initial = [-1, 0]
    xx = np.linspace(0, 1.57, 100)
    yarr = rk4(func = f, y0 = initial, a = 0, b = 1.57, n = 8)
    print(yarr)    
    # plt.plot(xx, result)
    # plt.show()

