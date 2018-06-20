#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 18:35:36 2018

@author: alex
----------------------------------
A property φ is transported by means of convection and diffusion through
the one-dimensional domain. The governing equation
is:
    
\frac{{\partial\phi}}{{\partial t}} + \frac{{\partial }}{{\partial x}}( \rho u  \phi) = \frac{{\partial }}{{\partial x}}\left(\Gamma\frac{{\partial  \phi}}{{\partial x}}\right)

The boundary conditions are:
    
                        φ(x, 0) = 0 para 0 ≤ x ≤ ∞
                        φ(0, t) = 1 para t > 0
                        φ(L, t) = 0 para t > 0, L → ∞
    
    
Compare the results with the analytical solution:

\phi \left ( x,t \right ) = \frac{1}{2}\left [ erfc\left ( \frac{x-ut}{2\sqrt{\Gamma t}} \right ) + exp\left ( \frac{ux}{D} \right )+erfc\left ( \frac{x+ut}{2\sqrt{\Gamma t}} \right )\right ]

___________________________________________________________

"""

#----------------------------------------------------------------------------------------
                               #Library Imports
#----------------------------------------------------------------------------------------

import FiniteVolumeMethod as fvm
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import  erfc
from matplotlib.animation import FuncAnimation

#----------------------------------------------------------------------------------------
                         #Problem Variables Declaration
#----------------------------------------------------------------------------------------

longitud = 2.5          # meters
phi_0 = 1               # intensive property at left boundary
phi_L = 0               # intensive property at right boundary
gamma = 0.001           # kg/ms
rho = 1.0                 # kg/m3
u = 1.0                   # m/s
N  = 350                 # Number of Nodes
dt = 0.002              # Time step (seconds)
tmax = 1                # Max time
iniCond = np.zeros(N+2) # Initial conditions of property at time zero

#----------------------------------------------------------------------------------------
                       #Mesh creation and mesh parameter definition
#----------------------------------------------------------------------------------------

malla = fvm.Mesh(nodes = N, length = longitud, timeInt = tmax, tsteps = int(tmax/dt))     # Mesh creation
nvx   = malla.volumes()                                                                   # Number of volumes
dx = malla.delta()                                                                        # Volumes size (grid step)
tSteps = malla.tsteps()                                                                   # Time steps

#----------------------------------------------------------------------------------------
                       #Problem variables printing
#----------------------------------------------------------------------------------------

fvm.printData(Longitud = longitud,
              Propiedad_izquierda = phi_0,
              Propiedad_derecha = phi_L,
              Conductividad = gamma,
              Densidad = rho,
              Velocidad = u,
              Nodos = N, 
              Volúmenes = nvx,
              Delta = dx, 
              Paso_de_Tiempo = dt,
              Tiempo_Máximo = tmax)
#----------------------------------------------------------------------------------------
                               #Object instances
#----------------------------------------------------------------------------------------

df1 = fvm.Diffusion1D(nvx, Gamma = gamma, dx = dx)     #Difussion Object
da1 = fvm.Advection1DCDS(nvx, rho = rho, dx = dx)    #Advection Object
dt1 = fvm.Temporal1DBackEuler(nvx,rho,dx,dt)           #Time Object

#----------------------------------------------------------------------------------------
                    #Time Advective-Diffusive coefficients allocation
#----------------------------------------------------------------------------------------
coef = fvm.Coefficients()       # Coefficients alias
coef.alloc(nvx)                 # Memory allocation

ux = np.zeros(nvx)              # Velocity vector allocation
ux += u                         # Velocity vectorization

#----------------------------------------------------------------------------------------
                    #Border conditions definition
#----------------------------------------------------------------------------------------

Sol = np.zeros(nvx)   # Initial array with zeroes
Sol[0]  = phi_0       # Left border condition (Dirichlet)
Sol[-1] = phi_L       # Right border condition (Dirichlet)

#----------------------------------------------------------------------------------------
                             #Spatial and Temporal grid allocation
#----------------------------------------------------------------------------------------

x = malla.createMesh()
t = malla.createTimeMesh()

#----------------------------------------------------------------------------------------
                             #Graph and animation initialization
#----------------------------------------------------------------------------------------
title_graf = '$\partial \phi / \partial t + \partial(p u \phi)/\partial x= \partial (\Gamma \partial\phi/\partial x)/\partial x$ con FVM'
fig = plt.figure(figsize=(8,4))                         # Figures
ax = plt.axes(xlim=(0, 3), ylim=(0, 1.2))               # Axis
line, = ax.plot(x, Sol, '--', label='FVM')
label = ax.text(2.6, 0.5, 'Time = {:>8.5f}'.format(0),
                ha='center', va='center',
                fontsize=12)

#----------------------------------------------------------------------------------------
                             #Analytycal Solution
#----------------------------------------------------------------------------------------
def analyticSol(x, u, t, Gamma):
 	divisor = 2 * np.sqrt(Gamma * t)
 	sol = 0.5 * (erfc((x - u * t)/ divisor) + 
 		np.exp(u * x) * np.exp(-Gamma) * erfc((x + u * t)/divisor))
 	return sol

#----------------------------------------------------------------------------------------
                            #Analytical solution graphication
#----------------------------------------------------------------------------------------

exac = analyticSol(x, u, dt * tSteps, gamma)
ax.plot(x,exac,'b-',label='Sol. Exac',lw=2)
plt.xlabel('$x$ [m]')
plt.ylabel('$\phi$ [...]')
ax.legend()

#----------------------------------------------------------------------------------------
                                #Function Solver
#----------------------------------------------------------------------------------------

def implicitSolver(i):
    time_step = i * dt
    #print('time step = {}'.format(time_step))

    coef.cleanCoefficients()                    # Initial values at 0

    df1.calcCoef()                              # Difussion coefficients calculation

    da1.setVel(ux)                              # Velocity definition

    da1.calcCoef()                              # Advection Coefficients
    
    da1.bcDirichlet('LEFT_WALL', Sol[0])        # Coefficients update (Left)
    da1.bcDirichlet('RIGHT_WALL', Sol[-1])      # Coefficients update (Right)

    dt1.calcCoef(Sol)                           # Temporal Coefficients
    
    #----------------------------------------------------------------------------------------
                        #Linear system solution matrix creation
    #----------------------------------------------------------------------------------------

    Su = coef.Su()                        # Source vector allocation
    A = fvm.Matrix(malla.volumes())       # Matrix allocation
    A.build(coef)                         # Matrix build
    
    #----------------------------------------------------------------------------------------
                                 #Linear system solution
    #----------------------------------------------------------------------------------------

    Sol[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])

    #----------------------------------------------------------------------------------------
                                 #Graph update
    #----------------------------------------------------------------------------------------
    line.set_ydata(Sol)                                                         # Y axis update
    label.set_text('Step = {:>8d} \n Time = {:>8.5f}'.format(i, time_step))
    ax.set_title(title_graf)

#----------------------------------------------------------------------------------------
                                #Solution animation
#----------------------------------------------------------------------------------------

anim = FuncAnimation(fig,               # La figura
                     implicitSolver,    # la función que cambia los datos
                     interval=1,        # Intervalo entre cuadros en milisegundos
                     frames=tSteps+1,   # Cuadros
                     repeat=False)       # Permite poner la animación en un ciclo

plt.show()
