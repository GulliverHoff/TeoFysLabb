# Python simulation of an electron in a 1d infinite box potential
# Integrate time independent SE using the Verlet method
# Locate eigenvalues by the shooting method
# MW 250402

import numpy as np
import matplotlib.pyplot as plt

h=6.62607015e-34   # Plancks constant
e=1.60217662e-19   # electron charge=-e

N=1000               # number of mesh points
dx=4/N             # step length
dx2=dx**2          # step length squared
c=2.0    # constant in Schrödinger equation

# potential energy function
def V(x):
    y = x**2/2 + x**4 # anharmonic oscillator
    return y



psi_initial, dpsi_initial = 1.0, 0.0 #even
# psi_initial, dpsi_initial = 0.0, 1.0 #odd

E_initial = 2.5
E_steps = 255
dE = 0.1

E_values = [E_initial + dE*i for i in range(E_steps)]
Psi_values = []
dPsi_values = []

for E in E_values:
    x = 0
    psi, dpsi = psi_initial, dpsi_initial

    x_tab = []          # list to store positions for plot
    psi_tab = []        # list to store wave function for plot
    dpsi_tab = []
    x_tab.append(x)
    psi_tab.append(psi)
    dpsi_tab.append(dpsi)

    for i in range(N) :
        d2psi = c*(V(x)-E)*psi
        psi += dpsi*dx + 0.5*d2psi*dx2
        d2psinew = c*(V(x+dx)-E)*psi
        dpsi += 0.5*(d2psi+d2psinew)*dx
        x += dx
        x_tab.append(x)
        psi_tab.append(psi)
        dpsi_tab.append(dpsi)
    Psi_values.append(psi)
    dPsi_values.append(dpsi)

zero_ranges = []


last_psi = 0
for i in range(len(E_values)):
    psi = Psi_values[i]
    if last_psi*psi < 0:
        zero_ranges.append(((E_values[i-1], Psi_values[i-1]), (E_values[i], Psi_values[i])))
    last_psi = psi

def subdivide(E_left, E_right, Ψ_left, Ψ_right, abs_Ψ_max, step=0) :
    E_mid = (E_left + E_right)/2
    x = 0
    Ψ = psi_initial
    dΨ = dpsi_initial
    for i in range(N):
        d2Ψ = c * (V(x) - E_mid) * Ψ
        Ψ += dΨ * dx + 0.5 * d2Ψ * dx2
        d2Ψnew = c * (V(x + dx) - E_mid) * Ψ
        dΨ += 0.5 * (d2Ψ + d2Ψnew) * dx
        x += dx
    Ψ_mid = Ψ
    if abs(Ψ_mid) < abs_Ψ_max or step > 900:
        return (E_mid, Ψ_mid)
    # print(f"{step} - E_left: {E_left} | E_right: {E_right} |---| Psi L:{Ψ_left} M:{Ψ_mid} R:{Ψ_right}")
    if Ψ_left * Ψ_mid < 0 : 
        return subdivide(E_left, E_mid, Ψ_left, Ψ_mid, abs_Ψ_max, step+1)
    else:
        return subdivide(E_mid, E_right, Ψ_mid, Ψ_right, abs_Ψ_max, step+1)

EigenValues = []
Psi_tolerance = 1e-5

for zero_range in zero_ranges:
    E_left = zero_range[0][0]
    E_right = zero_range[1][0]
    Psi_left = zero_range[0][1]
    Psi_right = zero_range[1][1]
    # print(f"E_left: {E_left}, E_right: {E_right}")

    EigenValues.append(subdivide(E_left, E_right, Psi_left, Psi_right, Psi_tolerance))

for pair in EigenValues:
    print(f"E: {pair[0]} | Psi: {pair[1]}")

# V_tab = []
# for x in x_tab:
#     V_tab.append(V(x))


# frame_y_max = max(max(psi_tab),max(dpsi_tab))
# frame_y_min = min(min(psi_tab),min(dpsi_tab))   

plt.close()
fig, ax = plt.subplots()
ax.plot(E_values, Psi_values, linewidth=2, color="#1B11E1")
ax.plot(E_values, dPsi_values, linewidth=2, color="#C71818")
ax.axhline(y=0, color='black', linewidth=1)
#plt.plot(x_tab, psi_tab, linewidth=2, linestyle="solid", color="#131FA0")
#plt.plot(x_tab, dpsi_tab, linewidth=2, linestyle="dotted", color="#A00008")
#plt.plot(x_tab, V_tab, linewidth=2, linestyle="dashed", color="#000000")
#plt.fill_between(x_tab, V_tab, color="#000000", alpha=0.3)

plt.ylim(-1e2, 1e2)
plt.xlabel('E',fontsize=15)
plt.ylabel('$\psi_E{(4)}$',fontsize=15)
#plt.savefig('psi.pdf')
plt.show()

