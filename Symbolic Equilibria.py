from sympy import roots, solve_poly_system
from sympy import *

# Define symbols
r, x, b, k, z, s_0, d_0, b_e, k_e, d_e, k_d, bet, y, v, g, a, c, q, sigma, gamm, p, h = symbols('r x b k z s_0 d_0 b_e k_e d_e k_d bet y v g a c q sigma gamm p h')

# Define parameter values as discussed in Section 5.1
b=1.02*10**(-9); a=1.333; c=1.8; q=100; d=2; gamm=1.5
sigma=1.83; p=2.4*10**(-4); h=5*10**4; g=10**5; s=5*10**3
bet=6*10**(-6); d_e=2; b_e=1, k_e=5*10**2; k_d=2.5*10**4
r =0.514

# Define model 2.1 in Wei paper
dxdt = r*x*(1-b*(x+y))*(g+x) - (bet*x*v)*(g+x)/(g+x) - k*x*z*(g+x)
dydt = (bet*x*v)*(g+x)/(g+x) - a*y*(g+x) - c*y*z*(g+x)
dvdt = q*a*y - sigma*v - gamm*v*z
dzdt = s_0*(h+y)*(k_d+x)*(k_e+x) - d_0*z*(h+y)*(k_d+x)*(k_e+x) + (b_e*x*z)*(h+y)*(k_d+x)*(k_e+x)/(k_e+x) - (d_e*x*z)*(h+y)*(k_d+x)*(k_e+x)/(k_d+x) + (p*y*z)*(h+y)*(k_d+x)*(k_e+x)/(h+y)

solve_poly_system([dxdt, dydt, dvdt, dzdt], x,y,v,z)