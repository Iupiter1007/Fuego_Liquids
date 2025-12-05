# mdot_vs_pressure.py
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq, newton

# -----------------------
# Fluid / geometry (changeable)
# -----------------------
rho = 786.0           # kg/m3 (IPA ~20°C)
mu = 0.00243          # Pa.s
Cd = 0.62             # discharge coeff for sharp orifice
D_or = 0.125 * 0.0254 # orifice diameter (1/8" -> m)
D_tube = 0.25 * 0.0254# tube ID (0.25" -> m)
L_tube = 0.3048       # 1 ft -> m
eps = 1e-6            # pipe absolute roughness (m), smooth nylon estimate
p_down = 101325.0     # plenum / downstream static pressure (Pa)

A_or = math.pi*(D_or/2)**2
A_tube = math.pi*(D_tube/2)**2

# -----------------------
# friction factor functions
# -----------------------
def haaland_f(Re, epsD):
    if Re <= 0:
        return 1.0
    if Re < 2300.0:
        return 64.0/Re
    inv_sqrt_f = -1.8 * math.log10((epsD/3.7)**1.11 + 6.9/Re)
    return 1.0/(inv_sqrt_f*inv_sqrt_f)

def colebrook_f_single(Re, epsD, tol=1e-12, maxiter=200):
    if Re < 2300.0:
        return 64.0/Re
    # initial guess: Haaland
    inv_sqrt = -1.8 * math.log10((epsD/3.7)**1.11 + 6.9/Re)
    f0 = 1.0/(inv_sqrt*inv_sqrt)
    # Newton on x = 1/sqrt(f)
    x = 1.0/math.sqrt(f0)
    for _ in range(maxiter):
        F = x + 2.0*math.log10(epsD/3.7 + 2.51*x/Re)
        dFdx = 1.0 + 2.0/math.log(10.0) * (2.51/Re) / (epsD/3.7 + 2.51*x/Re)
        dx = -F/dFdx
        x += dx
        if abs(dx) < tol:
            break
    return 1.0/(x*x)

# wrapper vectorized
def colebrook_f(Re, epsD):
    Re = np.asarray(Re, dtype=float)
    f = np.empty_like(Re)
    for i, Re_i in enumerate(Re):
        f[i] = colebrook_f_single(Re_i, epsD)
    return f

# -----------------------
# solver to get velocity in tube (v) for given upstream pressure and friction model
# Model: p_tank - p_down = 0.5*rho*(K_or + f*(L/D)) * v^2
# where K_or = (A_tube^2 / (A_or^2 * Cd^2))
# and v is tube velocity (Q/A_tube). Solve for v.
# -----------------------
def solve_v_for_pressure(p_tank, method='haaland'):
    # p_tank, p_down in Pa
    dP = p_tank - p_down
    if dP <= 0:
        return 0.0
    K_or = (A_tube**2)/(A_or**2) * (1.0/(Cd**2))
    epsD = eps / D_tube

    # residual function: 0.5*rho*(K_or + f(L/D))*v^2 - dP = 0
    def residual(v):
        if v <= 0:
            return -dP
        Re = rho * v * D_tube / mu
        if method == 'haaland':
            f = haaland_f(Re, epsD)
        elif method == 'colebrook':
            f = colebrook_f_single(Re, epsD)
        else:
            raise ValueError("method must be 'haaland' or 'colebrook'")
        return 0.5 * rho * (K_or + f * (L_tube / D_tube)) * v**2 - dP

    # bounds for v (m/s)
    v_low = 1e-8
    v_high = 1000.0  # very high upper bound; adjust if needed
    v = brentq(residual, v_low, v_high, xtol=1e-9, rtol=1e-9, maxiter=200)
    return v

# -----------------------
# sweep pressures and compute mdot
# -----------------------
def mdot_for_p(p_tank_pa, method='haaland'):
    v = solve_v_for_pressure(p_tank_pa, method=method)
    Q = A_tube * v
    return rho * Q, v

# default base pressure and doubled
p_base_psi = 800.0
p_base = p_base_psi * 6894.76
p_double = 2.0 * p_base

methods = ['haaland', 'colebrook']

results = {}
for method in methods:
    mdot_base, v_base = mdot_for_p(p_base, method=method)
    mdot_double, v_double = mdot_for_p(p_double, method=method)
    results[method] = {
        'p_base_pa': p_base,
        'p_double_pa': p_double,
        'mdot_base': mdot_base,
        'mdot_double': mdot_double,
        'v_base': v_base,
        'v_double': v_double,
        'ratio_mdot': mdot_double/mdot_base if mdot_base>0 else float('nan'),
        'percent_inc': (mdot_double/mdot_base - 1.0)*100.0 if mdot_base>0 else float('nan'),
    }

# Print summary
for method in methods:
    r = results[method]
    print(f"--- method: {method} ---")
    print(f"Base P = {p_base/6894.76:.1f} psi, Double P = {p_double/6894.76:.1f} psi")
    print(f"mdot_base = {r['mdot_base']:.6f} kg/s, mdot_double = {r['mdot_double']:.6f} kg/s")
    print(f"v_base = {r['v_base']:.3f} m/s, v_double = {r['v_double']:.3f} m/s")
    print(f"mdot ratio (double/base) = {r['ratio_mdot']:.6f} (~{r['percent_inc']:.2f}% increase)")
    print()

# -----------------------
# continuous curve and plot
# -----------------------
p_list_psi = np.linspace(100.0, 2000.0, 40)  # psi
p_list = p_list_psi * 6894.76

mdot_haal = []
mdot_cole = []
for p in p_list:
    mdot_haal.append(mdot_for_p(p, method='haaland')[0])
    mdot_cole.append(mdot_for_p(p, method='colebrook')[0])

plt.figure(figsize=(8,5))
plt.plot(p_list_psi, mdot_haal, '-o', markersize=4, label='Haaland')
plt.plot(p_list_psi, mdot_cole, '-s', markersize=4, label='Colebrook')
plt.axvline(p_base_psi, color='k', linestyle=':', label=f'base {p_base_psi:.0f} psi')
plt.axvline(2*p_base_psi, color='gray', linestyle='--', label=f'double {2*p_base_psi:.0f} psi')
plt.xlabel('Upstream tank pressure (psi)')
plt.ylabel('Mass flow (kg/s)')
plt.title('Mass flow vs upstream tank pressure (IPA, orifice+0.25\" tube, plenum at 1 atm)')
plt.grid(True, which='both', linestyle=':', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()
