# Python simulation of an electron in a 1d infinite box potential
# Integrate time independent SE using the Verlet method
# Locate eigenvalues by the shooting method
# MW 250402

import numpy as np
import matplotlib.pyplot as plt

e=1.60217662e-19   # electron charge=-e

N=1000               # number of mesh points
dx=4/N             # step length
dx2=dx**2          # step length squared
c=2.0   # constant in Schrödinger equation

# exact solution for infinite box potential
#print('E2=',EeV*2**2,'eV')

# input energy guess
#EeV = 0.3          # input energy in eV: test 0.3 , 0.4 , 0.3760 , 1.5
#E = EeV*e          # input energy in J

# potential energy function
def V(x):
    #y = 0.0
    #y = x**2/2 # harmonic oscillator
    y = x**2/2 + x**4 # anharmonic oscillator
    return y

# initial values and lists
x = 0               # initial value of position x

# even solution
psi = 1.0           # wave function at initial position
dpsi = 0.0          # derivative of wave function at initial position

# odd solution
# psi = 0.0           # wave function at initial position
# dpsi = 1.0          # derivative of wave function at initial position

x_tab = []          # list to store positions for plot
psi_tab = []        # list to store wave function for plot
x_tab.append(x)
psi_tab.append(psi)

#Odd
E = 2.7378797895563327
E = 7.942309864278653
E = 14.202840661135799
E = 21.23576792546728
#Even
E = 5.179250784485525


for i in range(N) :
    d2psi = c*(V(x)-E)*psi
    psi += dpsi*dx + 0.5*d2psi*dx2
    d2psinew = c*(V(x+dx)-E)*psi
    dpsi += 0.5*(d2psi+d2psinew)*dx
    x += dx
    x_tab.append(x)
    psi_tab.append(psi)

plt.close()
plt.plot(x_tab, psi_tab, linewidth=2)
plt.xlabel('x/a',fontsize=15)
plt.ylabel('$\psi$',fontsize=15)
#plt.savefig('psi.pdf')
plt.show()

