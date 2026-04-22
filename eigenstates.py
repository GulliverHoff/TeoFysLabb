# Python simulation of an electron in a 1d infinite box potential
# Integrate time independent SE using the Verlet method
# Locate eigenvalues by the shooting method
# MW 250402

import numpy as np
import matplotlib.pyplot as plt

h = 6.62607015e-34  # Plancks constant
hbar = h / (2 * np.pi)
m = 9.10938356e-31  # electron mass
e = 1.60217662e-19  # electron charge=-e

N = 10  # number of mesh points
a = 1.0e-9  # well width a=1 nm
dx = a / N  # step length
dx2 = dx**2  # step length squared
c = 2.0 * m / hbar**2  # constant in Schrödinger equation

# exact solution for infinite box potential
E = h**2 / (8 * m * a**2)  # Joule
EeV = E / e  # electron volt
print("E1=", EeV, "eV")
print("E2=", EeV * 2**2, "eV")

# input energy guess
# EeV = 0.3          # input energy in eV: test 0.3 , 0.4 , 0.3760 , 1.5
# E = EeV*e          # input energy in J


# potential energy function
def V(x):
    y = 0.0
    # y = x**2/2 # harmonic oscillator
    # y = x**2/2 + x**4 # anharmonic oscillator
    return y


# initial values and lists
x = 0  # initial value of position x

# even solution
# psi = 1.0           # wave function at initial position
# dpsi = 0.0          # derivative of wave function at initial position

# odd solution
psi = 0.0  # wave function at initial position
dpsi = 1.0  # derivative of wave function at initial position

x_tab = []  # list to store positions for plot
psi_tab = []  # list to store wave function for plot
x_tab.append(x / a)
psi_tab.append(psi)

for i in range(N):
    d2psi = c * (V(x) - E) * psi
    psi += dpsi * dx + 0.5 * d2psi * dx2
    d2psinew = c * (V(x + dx) - E) * psi
    dpsi += 0.5 * (d2psi + d2psinew) * dx
    x += dx
    x_tab.append(x / a)
    psi_tab.append(psi)

print("E=", EeV, "eV , psi(x=a)=", psi)

plt.close()
plt.plot(x_tab, psi_tab, linewidth=2)
plt.xlabel("x/a", fontsize=15)
plt.ylabel("$\psi$", fontsize=15)
# plt.savefig('psi.pdf')
plt.show()
