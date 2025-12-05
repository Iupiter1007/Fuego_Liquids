# requires scipy
import math
from scipy.optimize import brentq

# Inputs (change as needed)
p_tank_psi = 800.0
p_plenum = 101325.0       # set to atm if vented plenum
rho = 786.0               # kg/m3 for IPA ~20C
mu = 0.00243              # Pa.s
D_tube = 0.25*0.0254      # m (0.25")
L = 0.3048                # 1 ft in m
D_or = 0.125*0.0254       # m (orifice: 1/8")
Cd = 0.62

# conversions/areas
p_tank = p_tank_psi * 6894.76
A_tube = math.pi*(D_tube/2)**2
A_or = math.pi*(D_or/2)**2
K_or = (A_tube**2)/(A_or**2) * (1.0/(Cd**2))
eps = 1e-6

def haaland(Re, D):
    return (-1.8*math.log10((eps/D)/3.7 + 6.9/Re))**-2 if Re>0 else 1.0

def residual(v):
    Re = rho*v*D_tube/mu if v>0 else 1e-9
    f = haaland(Re, D_tube)
    return 0.5*rho*(K_or + f*L/D_tube)*v**2 - (p_tank - p_plenum)

v = brentq(residual, 1e-6, 500.0)
Re = rho*v*D_tube/mu
f = haaland(Re, D_tube)
Q = A_tube*v
mdot = rho*Q

print("v (m/s) =", v)
print("Q (m3/s) =", Q)
print("mdot (kg/s) =", mdot)
print("Re =", Re, "f =", f)
print("p_plenum (Pa) =", p_tank - 0.5*rho*(K_or + f*L/D_tube)*v**2)
