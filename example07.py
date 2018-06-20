#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 18:35:36 2018

@author: alex

Example 5.1 from Malalasekera Book
----------------------------------
A property φ is transported by means of convection and diffusion through
the one-dimensional domain sketched in Figure 5.2. The governing equation
is (5.3); the boundary conditions are φ 0 = 1 at x = 0 and φ L = 0 at x = L. Using
five equally spaced cells and the central differencing scheme for convection
and diffusion, calculate the distribution of φ as a function of x for (i) Case 1:
u = 0.2 m/s, (ii) Case 2: u = 2.5 m/s, and compare the results with the
analytical solution:

\frac{\phi -\phi_{0}}{\phi_{L}-\phi_{0}} = \frac{exp\left ( \rho ux/\Gamma  \right )-1}{exp\left ( \rho uL/\Gamma  \right )-1}

(iii) Case 3: recalculate the solution for u = 2.5 m/s with 20 grid nodes and
compare the results with the analytical solution. The following data apply:
length L = 1.0 m, ρ = 1.0 kg/m 3 , Γ = 0.1 kg/m.s.


___________________________________________________________

 |<-------------- L ----------------->|
 A                                    B
 |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
 |                                    |
 |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

phi = 1                           phi = 0

             -----------> u
___________________________________________________________

"""
#----------------------------------------------------------------------------------------
                               #Library Imports
#----------------------------------------------------------------------------------------

import FiniteVolumeMethod as fvm
import numpy as np
import matplotlib.pyplot as plt
from Coefficients import Coefficients

#----------------------------------------------------------------------------------------
                         #Problem Variables Declaration
#----------------------------------------------------------------------------------------

longitud = 1         # meters
phi_0 = 1            #intensive property at left boundary
phi_L = 0            #intensive property at right boundary
gamma = 0.1          # kg/ms
rho = 1              # kg/m3
u = 0.2              #m/s
N  = 6               # Number of nodes

#----------------------------------------------------------------------------------------
                       #Mesh creation and mesh parameter definition
#----------------------------------------------------------------------------------------

malla = fvm.Mesh(nodes = N, length = longitud)           # Mesh creation
nx    = malla.nodes()                                    # Number of nodes
nvx   = malla.volumes()                                  # Number of volumes
delta = malla.delta()                                    # Volumes size (grid step)

#----------------------------------------------------------------------------------------
                       #Problem variables printing
#----------------------------------------------------------------------------------------

fvm.printData(Longitud = longitud,
              Propiedad_izquierda = phi_0,
              Propiedad_derecha = phi_L,
              Conductividad = gamma,
              Densidad = rho,
              Velocidad = u,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)

#----------------------------------------------------------------------------------------
                               #Velocity allocation and vectorization
#----------------------------------------------------------------------------------------

ux = np.zeros(nvx-1)    # Domain velocity definition (used in costructor)
ux += u

#----------------------------------------------------------------------------------------
                               #Object instances
#----------------------------------------------------------------------------------------

df1 = fvm.Diffusion1D(nvx, Gamma = gamma, dx = delta)
da1 = fvm.Advection1DUpwind2(nvx, rho = rho, dx = delta, u = ux)

#----------------------------------------------------------------------------------------
                       #Advective-Diffusive coefficients allocation
#----------------------------------------------------------------------------------------

Coefficients.alloc(nvx)         # Memory allocation
df1.calcCoef()                  # Diffusice coefficients calculation
da1.calcCoef()                  # Advective coefficients calculation

#----------------------------------------------------------------------------------------
                    #Border conditions definition
#----------------------------------------------------------------------------------------

Sol = np.zeros(nvx)   # Initial array with zeroes
Sol[0]  = phi_0       # Left border condition (Dirichlet)
Sol[-1] = phi_L       # Right border condition (Dirichlet)

da1.bcDirichlet('LEFT_WALL', Sol[0])   # Coefficients update (Left)
da1.bcDirichlet('RIGHT_WALL', Sol[-1]) # Coefficients update (Right)

#----------------------------------------------------------------------------------------
                        #Linear system solution matrix creation
#----------------------------------------------------------------------------------------

Su = da1.Su()                       # Source vector allocation
A = fvm.Matrix(malla.volumes())     # Matrix allocation
A.build(df1)                        # Matrix build

#----------------------------------------------------------------------------------------
                                 #Linear system solution
#----------------------------------------------------------------------------------------

Sol[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
print('Solución = ', Sol)

#----------------------------------------------------------------------------------------
                             #Graph x axis initialization
#----------------------------------------------------------------------------------------

x = malla.createMesh()

#----------------------------------------------------------------------------------------
                             #Analytycal Solution
#----------------------------------------------------------------------------------------

phi_a = ((np.exp((rho*u*x)/gamma)-1)/(np.exp((rho*u*longitud)/gamma)-1))*(phi_L-phi_0) + phi_0

#----------------------------------------------------------------------------------------
                              #Results Plotting
#----------------------------------------------------------------------------------------

x *= 100                                                                        # [cm] transformation
plt.plot(x,phi_a, '-', label = 'Sol. analítica')                                # Analytical Solution curve
plt.plot(x,Sol,'o', label = 'Sol. FVM')
plt.title('Ecuación de Advección Difusión Upwind segundo orden [example07]')
plt.xlabel('$x$ [cm]')
plt.ylabel('Propiedad (phi) ')
plt.grid()
plt.legend()
plt.savefig('example07.pdf')
plt.show()
