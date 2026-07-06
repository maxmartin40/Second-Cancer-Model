from sympy import *
import numpy as np
from sympy import init_printing

# Define symbols
# Beta = bet
# alpha = a
# gamma = gamm
# sigma = d
r, x, b, k, z, s_0, d_0, b_e, k_e, d_e, k_d, bet, y, v, g, a, c, q, d, gamm, p, h = symbols('r x b k z s_0 d_0 b_e k_e d_e k_d bet y v g a c q d gamm p h')

# Define variables that have been calcuated by hand
x_1 = (r-b*r*y-bet**2*v-c*k*z)
y_1 = (-bet**2*r+bet**4*v+c*k*bet**2*z+c**2*r*z-c**2*bet**2*v*z-c**3*k*z**2)/(b*r*(-bet**2-a+c**2*z))
numerator_v = a*bet**2*r*q-bet**2*r-c*a*k*bet**2*z*q+c*k*bet**2*z-c**2*a*r*z*q+c**2*r*z+c**3*a*k*z**2*q-c**3*k*z**2-c*q*bet**2*r*z+c**2*g*k*bet**2*z**2-c*a*g*r*z+c**2*a*g*k*z**2+c**3*g*r*z**2-c**4*g*k*z**3
denominator_v = a*bet**4*q-bet**4-c**2*a*bet*z*q+c**2*bet*z+b*d*bet**2*r+a*b*d*r-c**2*b*d*r*z-c*gamm*bet**4*z-c*a*gamm*bet**2*z+c**3*gamm*bet**2*z**2
v_1 = numerator_v/denominator_v
dz_dt = s_0 - d_0*z + (b_e*x*z)/(k_e+x) - (d_e*x*z)/(k_d+x) + (p*y*z)/(h+y)

# Substitute in new expressions for x, y, and v
dz_dt_a = dz_dt.subs(x,x_1)
dz_dt_b = dz_dt_a.subs(y,y_1)
dz_dt_c = dz_dt_b.subs(v,v_1)

# dz_dt_c now represents Equation 4 (see work) that is fully substituted in so we need to try and factor and get a polynomial

print(simplify(expand(dz_dt_c)))