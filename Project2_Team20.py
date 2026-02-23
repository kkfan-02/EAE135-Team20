# Project 2 - EAE 135
# Team 20
# DUE: 02/24/26

import numpy as np
import matplotlib.pyplot as plt

F = 726000          # max thrust of first stage (N)
rocket_diam = 1.28  # rocket diameter (m)
prop_diam = 1.18    # propellant diameter (m)
AoA = 20            # angle of attack (deg)
length = 5*rocket_diam      # rocket length (m)
thickness_oneLayer = 0.050 / 8        # thickness of one layer of carbon epoxy (m)
cpcg = rocket_diam  # distance btwn cp-cg (m)
P1_net = 600000          # net axial force on O (N)
Weight = 8500            # weight of rocket (kg)
lift = 1400000      # lift force (N)
FoS = 1.25          # spacecraft safety factor

# AS4/Epoxy properties
rho = 1522.3948             # density (kg/m^3)
StrengthTensile_0deg = 2137  # MPa - Axial Tensile Strength (0 deg ply)
StrengthCompressive_0deg = -184 * 6.89476 # -184 ksi to MPa (0 deg ply)
StrengthTensile_90deg = 53.4 # Mpa
StrengthCompressive_90deg = -24.4 * 6.89476 # -24.4 ksi to MPa (90 deg ply)
StrengthTensile_45deg = 53.4 # Mpa
StrengthCompressive_45deg = -24.4 * 6.89476 # -24.4 ksi to MPa (+/-45 deg ply)

# Propellant Properties (HTPB, Aged 185 Days)



"""
# Axial Stress/Displacement
A = np.pi*(rocket_diam/2)**2 - np.pi*(prop_diam/2)**2  # cross-sectional area of the case (m^2)

E0  = E11
E90 = E22
E45 = 1/(1/E11 + 1/E22)
Eeff = (2*E0 + 4*E45 + 2*E90)/8        # Simple effective modulus for [0/+45/-45/90]s *idk maybe*

P_tot = P - lift*np.sin(np.radians(AoA))
u1_axial = lambda x1: (P_tot/(Eeff*A))*(x1 - length)

u_left = u1_axial(0)           # max magnitude
u_right = u1_axial(length)     # fixed end

x1 = np.linspace(0, length, 200)
u_axial = u1_axial(x1)

plt.figure()
plt.plot(x1, u_axial, linewidth=1.5)
plt.grid(True)
plt.xlabel('x_1 (m)')
plt.ylabel('u_{1,axial} (m)')
plt.title('Axial displacement along beam')
plt.show()
"""
# 1.) Find Bending Stiffness H33_C
# Layup [0/+/-45/90]s
# ----------------- [0 deg]
# ----------------- [+45 deg]
# ----------------- [-45 deg]
# ----------------- [90 deg]
# ----------------- [90 deg]
# ----------------- [-45 deg]
# ----------------- [+45 deg]
# ----------------- [0 deg]
# \\\\\\\\\\\\\\\\\ [CORE]




# 2.) Iterating Layer Integrals


# 3.) Compute Q Matrix


# 4a.) Find Qbar (Iterate over orientation angle)


# 4b.) Sbar

# 4c.) Extract Projected Young's Modulus for each layer


# 5.) 