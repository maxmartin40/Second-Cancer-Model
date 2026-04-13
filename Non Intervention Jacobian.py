from sympy import *
import numpy as np

# Define symbols
r, x, b, k, z, s_0, d_0, b_e, k_e, d_e, k_d, lamb = symbols('r x b k z s_0 d_0 b_e k_e d_e k_d lamb')

# Define system
dx_dt = r*x*(1-b*x)-k*x*z
dz_dt = s_0 - d_0*z + (b_e*x*z)/(k_e+x) - (d_e*x*z)/(k_d+x)

# First row of Jacobian
dx_dt_dx = diff(dx_dt, x)
dx_dt_dz = diff(dx_dt, z)

# Second row of Jacobian
dz_dt_dx = diff(dz_dt, x)
dz_dt_dz = diff(dz_dt, z)

#print(simplify(dx_dt_dx))

# Parameter values
r=0.514; b=1.02*10**(-9); s_0=5000; k=10**(-4)
#d_0=4.1; b_e=5.5; d_e=5.50001; k_e=1100; k_d=24000 #A and B in Figure 1
d_0=5; b_e=6; d_e=4; k_e=1100; k_d=24000 #C and D in Figure 2

def J(x,z):
    A_11 = -b*r*x - k*z + r*(-b*x + 1)
    A_12 = -k*x
    A_21 = -(b_e*x*z)/(k_e + x)**2 + (b_e*z)/(k_e + x) + (d_e*x*z)/(k_d + x)**2 - d_e*z/(k_d + x)
    A_22 = (b_e*x)/(k_e + x) - d_0 - (d_e*x)/(k_d + x)
    A = Matrix([[A_11-lamb,A_12],[A_21,A_22-lamb]])
    return det(A)

'''
Equilibria for A and B from Figure 1
Ensure the correct parameters are commented out in Lines 23-24
'''
#x_0, z_0 = (2004.54, 5139.99) #E_1
#x_0, z_0 = (13170.77, 5139.93) #E_2
#x_0, z_0 = (747776127.28, 1219.56) #E_3
#x_0, z_0 = (0,1219.51) #E_4

#print(J(x_0,z_0))

# Coefficient numbers match equation numbers in Overleaf file
# coefficients_7 = [1,  0.972724064268509, 0.445634059727609]
# coefficients_8 = [1, 0.972776396941385, -0.445640836780393]
# coefficients_9 = [1, 4.4918857085892, 1.60729830709304]
# coefficients_10 = [1, 3.707951, -1.6074009]

# print(np.roots(coefficients_4))

'''
Equilibria for C and D from Figure 1
Ensure the correct parameters are commented out in Lines 23-24
'''
#x_0, z_0 = (3351.47, 5139.98) #E_1
#x_0, z_0 = (15648.35, 5139.92) #E_2
#x_0, z_0 = (662481567.52, 1666.74) #E_3
#x_0, z_0 = (0, 1000) #E_4

#print(J(x_0,z_0))

# coefficient numbers match equations in Overleaf file
coefficients_11 = [1, 0.972791607965718, 0.352707926013966]
coefficients_12 = [1, 0.972790945979541, -0.30193257399118]
coefficients_13 = [1, 3.34719073052093, 1.04190765808692]
coefficients_14 = [1, 4.586, -2.07]

#print(np.roots(coefficients_14))