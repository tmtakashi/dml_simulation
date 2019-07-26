'''
Implementation of
"Flat-Panel Loudspeaker Simulation Model with Electromagnetic Inertial Exciters and Enclosures"
by David A. Anderson
https://doi.org/10.17743/jaes.2017.0027
'''
import numpy as np
from util import mode_function

'''
Parameters
'''

R = 25 # number of in each dimention
R_tot = R ** 2 # total number of modes




# -------- Panel ----------
E = 3.2e+9 # Young's modulus of the panel
h = 3e-3 # thickness of the panel
nu = 0.5 # panel material Poisson's ratio
rho = 1180 # density of the panel
Ip = h**3 / (12 * (1 - nu**2)) # panel moment of inertia 
D = E * h**3 / (12*(1 - nu**2)) # panel bending stiffness

# dimensions
Lx = 0.35
Ly = 0.245

rm = np.arange(1, R + 1) # mode index in x dimention
rn = np.arange(1, R + 1) # mode index in y dimention

deltam = 1 / ((rn * Lx / (rm * Ly))**2 + 2) # horizontal edge effect factor
deltan = 1 / ((rm * Ly / (rn * Lx))**2 + 2) # vertical edge effect factor

Qr = 9 # quality factor of each mode
br = 0.5 / Qr # damping ratio of each mode

beta = Lx * Ly / 4 # modal effective pistonic area

Phir = (np.pi**2 * rm**2 / Lx**2) + (np.pi**2 * rn**2 / Ly**2)

Zpr = complex(D * Phir**2, omega * br - omega**2 * rho * h) # panel mode r impedance wrt displacement

# ------- Exiters --------
I = 1 # total number of exciter
mi = np.ones(I) * 0.032 # mass of each exiter
ki = np.ones(I) * 2000 # spring constant of each exiter
Qi = np.ones(I) * 9 # quality factor of each exciter
bi = 0.5 / Qi # damping ratio of each exiter
Ftildai = np.ones(I) * 0.01 # output force of each exciter
omega = np.sqrt(ki / mi) # resonant angular frequency of each exiter

# exciter compliance and resistance wrt displacement
Zci = np.array([complex(ki[i], omega[i] * bi[i]) for i in range(I)])
# exciter compliance and resistance wrt displacement
Zmi = np.array([complex(ki[i], omega[i] * bi[i] - omega[i]**2 * mi[i]) for i in range(I)])

xi = 0.37 * Lx
yi = 0.69 * Ly

alphair = np.zeros((R_tot, I))
for r in range(R_tot):
    for i in range(I):
        


# ------ Enclosure -------
c0 = 340 # sound speed of air
rho0 = 1.29 # density of air
V = Lx * Ly * 0.0254


Cr = np.zeros((R, R)) # acoustic coupling factor
for i in range(R):
    for j in range(R):
        if i * j % 2 != 0:
            Cr[i, j] = 4 / (i * j * np.pi**2)

kb = (Lx * Ly)**2 * rho0 * c0**2 / V

# ----- Panel Mode -------
W = np.zeros((R, R)) # panel mode
for m in range(R):
    for n in ramge(R):
        W[m, n] = mode_function(xi, yi, m, n, Lx, Ly)

omegar = np.zeros((R, R)) # resonant frequency of each mode
for i in range(R):
    for j in range(R):
        phitilda = ((i + deltam[i]) / Lx)**2 + ((i + deltan[j]) / Ly)**2
        omegar[i, j] = np.pi**2 * np.sqrt(D / (rho * h)) * phitilda

order_indices = omegar.ravel().argsort() # get indices in ascending order

W_flat = W.ravel()[order_incedices] # get mode in order of resonant frequency

'''
Create matrix B
'''
B = np.zeros((R, R))
for i in range(R):
    for j in range(R):
        if i == j:
            B[i, j] = beta * Zpr[i] + Cr[i, j]**2 * kb + 
