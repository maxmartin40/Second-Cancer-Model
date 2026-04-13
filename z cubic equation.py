import numpy as np
from sympy import *
import sympy as sp
import pprint
import numpy.polynomial.polynomial as poly

# Define symbols
r, x, b, k, z, s_0, d_0, b_e, k_e, d_e, k_d, lamb = symbols('r x b k z s_0 d_0 b_e k_e d_e k_d lamb')

# x = (r-k*z)/(b*r)
dz_dt = s_0 - d_0*z + (b_e*(r-k*z)/(b*r)*z)/(k_e+(r-k*z)/(b*r)) - (d_e*(r-k*z)/(b*r)*z)/(k_d+(r-k*z)/(b*r))

# simp = sp.simplify(dz_dt)
# common_denom = sp.together(simp)
# num, denom = sp.fraction(common_denom)
# num_expand = sp.expand(num)
# factored = sp.factor(num_expand)
# pprint.pprint(factored)

# Parameter values
r=0.514; b=1.02*10**(-9); s_0=5000; k=10**(-4)
#d_0=4.1; b_e=5.5; d_e=5.50001; k_e=1100; k_d=24000 #A and B in Figure 1
d_0=5; b_e=6; d_e=4; k_e=1100; k_d=24000 #C and D in Figure 2

'''
This uses the work that SymPy did in lines 13-18 to factor the z equation
'''
# number of coefficient variable is the power of z so for example: coeff0 is the z^0 term
coeff0 = r**2*s_0 + b*k_d*r**2*s_0 + b*k_e*r**2*s_0 + b**2*k_d*k_e*r**2*s_0 #c
coeff1 = b_e*r**2 - d_0*r**2 - d_e*r**2 + b*b_e*k_d*r**2 - b*d_0*k_d*r**2 - b*d_0*k_e*r**2 - b*d_e*k_e*r**2 - b**2*d_0*k_d*k_e*r**2 - 2*k*r*s_0 - b*k*k_d*r*s_0 - b*k*k_e*r*s_0 #z^1
coeff2 = -2*b_e*k*r + 2*d_0*k*r + 2*d_e*k*r - b*b_e*k*k_d*r + b*d_0*k*k_d*r + b*d_0*k*k_e*r + b*d_e*k*k_e*r + k**2*s_0 #z^2
coeff3 = b_e*k**2 - d_0*k**2 - d_e*k**2 #z^3

# Two different ways to run the roots equation and then printing them, notice that they take them in reverse order
z_2_coefficients = [coeff3, coeff2, coeff1, coeff0]
z_2_coefficients_poly = [coeff0, coeff1, coeff2, coeff3]
z_1, z_2, z_3 = np.roots(z_2_coefficients)
print([z_1,z_2,z_3])
#print(np.roots(z_2_coefficients))
#print(z_2_coefficients)

# Solve for corresponding x values
x_1 = (r-k*z_1)/(b*r)
x_2 = (r-k*z_2)/(b*r)
x_3 = (r-k*z_3)/(b*r)

print([x_1,x_2,x_3])