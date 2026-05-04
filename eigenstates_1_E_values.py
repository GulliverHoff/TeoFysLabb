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

N = 1000  # number of mesh points
a = 1.0e-9  # well width a=1 nm
dx = a / N  # step length
dx2 = dx**2  # step length squared
c = 2.0 * m / ħ**2  # constant in Schrödinger equation

for n in range(1, 6):
    # exact solution for infinite box potential
    E_n = (h * n / a)**2 / (8 * m)  # Joule
    EeV = E_n / e  # electron volt
    print(f"E{sub(n)}=", EeV, "eV")
#print(f"E{sub(2)}=", EeV * 2**2, "eV")

# input energy guess
EeV_in = 10  # input energy in eV: test 0.3 , 0.4 , 0.3760 , 1.5
E_in = EeV_in * e  # input energy in J
dE = 0.1 * e

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


# def subdivide(x_low, x_high, Ψ_left, dx_threshold):
#     dx = x_high - x_low
#     d2Ψ = c * (V(x) - E_in) * Ψ_left
#     Ψ = Ψ_left + dΨ * dx / 2 + 0.5 * d2Ψ * dx2 / 4

    # if last_Ψ * Ψ < 0 :
    #     zero_regions.append([x, x+dx])


# last_Ψ = 0

E = 0
Ψa_last = 0
E_last = E
EΨ_zero_spans = []
while E < E_in:
    E += dE
    Ψ = 0.0  # wave function at initial position
    dΨ = 1.0  # derivative of wave function at initial position

    x_tab = []  # list to store positions for plot
    Ψ_tab = []  # list to store wave function for plot
    x_tab.append(x / a)
    Ψ_tab.append(Ψ)

    for i in range(N):
        d2Ψ = c * (V(x) - E) * Ψ
        last_Ψ = Ψ
        Ψ += dΨ * dx + 0.5 * d2Ψ * dx2
        d2Ψnew = c * (V(x + dx) - E) * Ψ
        dΨ += 0.5 * (d2Ψ + d2Ψnew) * dx
        x += dx
        x_tab.append(x / a)
        Ψ_tab.append(Ψ)
    if Ψa_last * Ψ < 0:
        EΨ_zero_spans.append([(E_last, Ψa_last), (E, Ψ)])
    Ψa_last = Ψ
    E_last = E

# for span in EΨ_zero_spans:
#     print(f"E_left: {round(span[0][0]/e, 2)} | E_right: {round(span[1][0]/e, 2)} - Ψ_left: {span[0][1]} | Ψ_right: {span[1][1]}")

def subdivide(E_left, E_right, Ψ_left, Ψ_right, subdiv_step) :
    E_mid = (E_left + E_right)/2
    x = 0
    Ψ = 0
    dΨ = 1.0
    for i in range(N):
        d2Ψ = c * (V(x) - E_mid) * Ψ
        Ψ += dΨ * dx + 0.5 * d2Ψ * dx2
        d2Ψnew = c * (V(x + dx) - E_mid) * Ψ
        dΨ += 0.5 * (d2Ψ + d2Ψnew) * dx
        x += dx
    Ψ_mid = Ψ
    if subdiv_step < 1:
        return (E_mid, Ψ_mid)
    if Ψ_left * Ψ_mid < 0 : 
        return subdivide(E_left, E_mid, Ψ_left, Ψ_mid, subdiv_step - 1)
    else:
        return subdivide(E_mid, E_right, Ψ_mid, Ψ_right, subdiv_step - 1)


EΨ_values = []
for i in range(len(EΨ_zero_spans)):
    span = EΨ_zero_spans[i]
    E_left = span[0][0]
    E_right = span[1][0]
    Ψ_left = span[0][1]
    Ψ_right = span[1][1]
    EΨ_values.append(subdivide(E_left, E_right, Ψ_left, Ψ_right, 40))
    
    
for EΨ in EΨ_values:
    print(f"E: {EΨ[0]/e} | Ψ: {EΨ[1]}")





# print(f"E{sub(n)}= {EeV_in} eV\n\u03A8(a)= {Ψ}")

# plt.close()
# fig, ax = plt.subplots()
# ax.plot(x_tab, Ψ_tab, linewidth=2)
# ax.axhline(y=0, color='black', linewidth=1)
# for x_start, x_end in zero_regions:
#     ax.fill_between([x_start/a, x_end/a], 0, 1,
#                     color='green', alpha=0.5, transform=ax.get_xaxis_transform())
# plt.xlabel("x/a", fontsize=15)
# plt.ylabel("$\psi$", fontsize=15)
# # plt.savefig('psi.pdf')
# plt.show()



