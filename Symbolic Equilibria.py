from sympy import roots, solve_poly_system
from sympy import *

# Define symbols
r, x, b, k, z, s_0, d_0, b_e, k_e, d_e, k_d, B, y, v, g, a, c, q, sigma, gamm, p, h = symbols('r x b k z s_0 d_0 b_e k_e d_e k_d B y v g a c q sigma gamm p h')

# Define model 2.1 in Wei paper
dxdt = r*x*(1-b*(x+y)) - (B*x*v)/(g+x) - k*x*z
dydt = (B*x*v)/(g+x) - a*y - c*y*z
dvdt = q*a*y - sigma*v - gamm*v*z
dzdt = s_0 - d_0*z + (b_e*x*z)/(k_e+x) - (d_e*x*z)/(k_d+x) + (p*y*z)/(h+y)