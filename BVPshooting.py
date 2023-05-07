import numpy as np
import matplotlib.pyplot as plt

def function(X):
	derivative = np.zeros(2)
	derivative[0] = X[1]
	derivative[1] = 2*(X[0])**3
	return derivative

def rk4(IC,x_array,N):
	#Step Size
	h = (x_array[-1] - x_array[0])/N
	x_o = x_array[0] ; y_o = IC[0] ; z_o = IC[1]
	y_array = [y_o] ; z_array = [z_o]
	
	for i in range(N):
		X = [y_o,z_o,x_o]
		k1 = h * function(X)
		X = [y_o + k1[0]/2,z_o + k1[1]/2,x_o + h/2]
		k2 = h * function(X)
		X = [y_o + k2[0]/2,z_o + k2[1]/2,x_o + h/2]
		k3 = h * function(X)
		X = [y_o + k3[0],z_o + k3[1],x_o + h]
		k4 = h * function(X)
		
		y_o = y_o + 1/6 * (k1[0] + 2*(k2[0] + k3[0]) + k4[0])
		z_o = z_o + 1/6 * (k1[1] + 2*(k2[1] + k3[1]) + k4[1])
		
		y_array.append(y_o) ; z_array.append(z_o)
		x_o = x_o + h
		
	return y_array,z_array
	
def secant_method(s0,s1,phi): # phi = [phi(s_i),phi(s_i+1)]
	s2 = s1 - ((s1-s0)/(phi[1]-phi[0]))*phi[1]
	return s2
	
def shooting_method(x_array,s0,s1,N,N_max,tol,eq):

	if eq == 1:
		phi = np.zeros(2)

		IC_1 = [1/3,s0] ; IC_2 = [1/3,s1]
		
		y_s0_arr,z_s0_arr = rk4(IC_1,x_array,N)
		y_s0 = y_s0_arr[-1]
		
		y_s1_arr,z_s1_arr = rk4(IC_2,x_array,N)
		y_s1 = y_s1_arr[-1]
		
		phi[0] = (y_s0 - 1/4)
		phi[1] = (y_s1 - 1/4)
		
		for i in range(N_max):
			if abs(phi[0]) <= tol:
				break
			else:
				s2 = secant_method(s0,s1,phi)
				
				s0 = s1
				s1 = s2
				IC_3 = [1/3,s1]
				
				y_s2_arr,z_s2_arr = rk4(IC_3,x_array,N)
				y_s2 = y_s2_arr[-1]
				
				phi[0] = phi[1]
				phi[1] = (y_s2 - 1/4)
	elif eq == 2:
		phi = np.zeros(2)

		IC_1 = [s0,-1/9] ; IC_2 = [s1,-1/9]
		
		y_s0_arr,z_s0_arr = rk4(IC_1,x_array,N)
		z_s0 = z_s0_arr[-1]
		
		y_s1_arr,z_s1_arr = rk4(IC_2,x_array,N)
		z_s1 = z_s1_arr[-1]
		
		phi[0] = (z_s0 - (-1/16))
		phi[1] = (z_s1 - (-1/16))
		
		for i in range(N_max):
			if abs(phi[0]) <= tol:
				break
			else:
				s2 = secant_method(s0,s1,phi)
				
				s0 = s1
				s1 = s2
				IC_3 = [s1,-1/9]
				
				y_s2_arr,z_s2_arr = rk4(IC_3,x_array,N)
				z_s2 = z_s2_arr[-1]
				
				phi[0] = phi[1]
				phi[1] = (z_s2 - (-1/16))
				
	elif eq == 3:
		phi = np.zeros(2)

		IC_1 = [(2-(9*s0))/3,s0] ; IC_2 = [(2-(9*s1))/3,s1]
		
		y_s0_arr,z_s0_arr = rk4(IC_1,x_array,N)
		y_s0 = y_s0_arr[-1]
		
		y_s1_arr,z_s1_arr = rk4(IC_2,x_array,N)
		y_s1 = y_s1_arr[-1]
		
		phi[0] = (y_s0 - 1/4)
		phi[1] = (y_s1 - 1/4)
		
		for i in range(N_max):
			if abs(phi[0]) <= tol:
				break
			else:
				s2 = secant_method(s0,s1,phi)
				
				s0 = s1
				s1 = s2
				IC_3 = [(2-(9*s1))/3,s1]
				
				y_s2_arr,z_s2_arr = rk4(IC_3,x_array,N)
				y_s2 = y_s2_arr[-1]
				
				phi[0] = phi[1]
				phi[1] = (y_s2 - 1/4)
			
	return y_s2_arr,z_s2_arr
	
def analytic(x_array):
	sol_anal = 1/(x_array + 3)
	return sol_anal
	
if __name__ == "__main__":
	N = 100
	x_array = np.linspace(0,1,N+1)
	s0,s1 = input("Enter the two guess (seperate them with ,): ").split(",")
	s0 = float(s0) ; s1 = float(s1)
	
	y,z = shooting_method(x_array,s0,s1,N,100,10**(-2),3)
	sol_anal = analytic(x_array)
	print(sol_anal)
	
	plt.plot(x_array,y)
	plt.plot(x_array,sol_anal,ls = "--")
	plt.show()

	

