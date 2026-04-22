# Python simulation of an electron in a 1d infinite box potential
# Integrate time independent SE using the Verlet method
# Locate eigenvalues by the shooting method
# MW 250402
from print_util import sub, sup
import numpy as np
import matplotlib.pyplot as plt

h = 6.62607015e-34  # Plancks constant
ħ = h / (2 * np.pi)
m = 9.10938356e-31  # electron mass
e = 1.60217662e-19  # electron charge=-e

N = 10  # number of mesh points
a = 1.0e-9  # well width a=1 nm
dx = a / N  # step length
dx2 = dx**2  # step length squared
c = 2.0 * m / ħ**2  # constant in Schrödinger equation

n = 1
# exact solution for infinite box potential
E_n = (h * n / a)**2 / (8 * m)  # Joule
EeV = E_n / e  # electron volt
print(f"E{sub(1)}=", EeV, "eV")
#print(f"E{sub(2)}=", EeV * 2**2, "eV")

# input energy guess
EeV_in = 1.5  # input energy in eV: test 0.3 , 0.4 , 0.3760 , 1.5
E_in = EeV_in * e  # input energy in J

# potential energy function
def V(x):
    y = 0.0
    # y = x**2/2 # harmonic oscillator
    # y = x**2/2 + x**4 # anharmonic oscillator
    return y


# initial values and lists
x = 0  # initial value of position x

# even solution
# Ψ = 1.0           # wave function at initial position
# dΨ = 0.0          # derivative of wave function at initial position

# odd solution
Ψ = 0.0  # wave function at initial position
dΨ = 1.0  # derivative of wave function at initial position

x_tab = []  # list to store positions for plot
Ψ_tab = []  # list to store wave function for plot
x_tab.append(x / a)
Ψ_tab.append(Ψ)

for i in range(N):
    d2Ψ = c * (V(x) - E_in) * Ψ
    Ψ += dΨ * dx + 0.5 * d2Ψ * dx2
    d2Ψnew = c * (V(x + dx) - E_in) * Ψ
    dΨ += 0.5 * (d2Ψ + d2Ψnew) * dx
    x += dx
    x_tab.append(x / a)
    Ψ_tab.append(Ψ)

print(f"E{sub(n)}= {EeV_in} eV\n\u03A8(a)= {Ψ}")

plt.close()
plt.plot(x_tab, Ψ_tab, linewidth=2)
plt.xlabel("x/a", fontsize=15)
plt.ylabel("$\psi$", fontsize=15)
# plt.savefig('psi.pdf')
#plt.show()
