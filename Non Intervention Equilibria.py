from sympy import *
import numpy as np
import random
from scipy.optimize import root
from scipy.optimize import fsolve, minimize
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Parameter values
r=0.514; b=1.02*10**(-9); s_0=5000; k=10**(-4)
#d_0=4.1; b_e=5.5; d_e=5.50001; k_e=1100; k_d=24000 #A and B in Figure 1
d_0=5; b_e=6; d_e=4; k_e=1100; k_d=24000 #C and D in Figure 2

def sys(t, z):
    x, y = z
    # dx = r*x*(1-b*x)-k*x*z
    # dz = s_0 - d_0*z + (b_e*x*z)/(k_e+x) - (d_e*x*z)/(k_d+x)
    return [r*x*(1-b*x)-k*x*z, s_0 - d_0*z + (b_e*x*z)/(k_e+x) - (d_e*x*z)/(k_d+x)]

sol = solve_ivp(sys, [0,15], [10,5], dense_output=True)

t = np.linspace(0, 15, 300)
z = sol.sol(t)
plt.plot(t, z.T)
plt.xlabel('t')
plt.legend(['x','z'], shadow=True)
plt.title('Lotka-Volterra System')
plt.show()