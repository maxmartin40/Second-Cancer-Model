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

def F(z, k, d_e):
    qty = s_0 - d_0*z + (b_e*(r-k*z)/(b*r)*z)/(k_e+(r-k*z)/(b*r)) - (d_e*(r-k*z)/(b*r)*z)/(k_d+(r-k*z)/(b*r))
    return qty

def find_equilibria_2(num_starts=25000, bounda=1000, boundb=5500):
    # Initialize lists
    roots = []
    message_list = []
    d_e_y = []
    k_x = []
    # Generate random guesses inside range from above
    for i in range(num_starts):
        k = random.uniform(10**(-5), 10**(-2))
        d_e = random.uniform(10**(-2), 1000)
        guess = random.uniform(bounda,boundb)
        # Run a solution to the system
        sol = root(F, x0=guess, method='lm',options={'xtol': 1e-14,'ftol': 1e-14}, args=(k, d_e))
        # Determine if there is a duplicate root and if not, add to root list
        if sol.success==True:
            r = sol.x
            solution = round(r[0])
            if solution not in roots:
                roots.append(r)
                k_x.append(k)
                d_e_y.append(d_e)
    #print(roots)
    return (k_x, d_e_y, message_list)

equilibria = find_equilibria_2()

# Plot solutions
d_e_y = equilibria[1]
k_x = equilibria[0]
plt.plot(k_x,d_e_y, '.')
plt.xlabel("k")
plt.ylabel("d_e")
plt.title("Recreation of Figure A1 (Garcia paper)")
plt.grid(True)
plt.show()