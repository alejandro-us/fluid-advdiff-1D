#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 15:32:38 2018

@author: luiggi

Example 4.2 from Malalasekera Book
----------------------------------
Now we discuss a problem that includes sources other than those arising
from boundary conditions. Figure 4.6 shows a large plate of thickness
L = 2 cm with constant thermal conductivity k = 0.5 W/m.K and uniform
heat generation q = 1000 kW/m3. The faces A and B are at temperatures
of 100°C and 200°C respectively. Assuming that the dimensions in the y- and
z-directions are so large that temperature gradients are significant in the 
x-direction only, calculate the steady state temperature distribution. Compare
the numerical result with the analytical solution. The governing equation is

$
\dfrac{d}{dx} \left( k \dfrac{d T}{dx} \right) + q = 0
$

___________________________________________________________

 |<-------------- 2.0 m ------------->|
 A                                    B
 |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
 |                                    |
 |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

T_A = 100                            T_B = 200
___________________________________________________________
Figure 4.6

Analytical solution:
    
$T = \left( \frac{T_B - T_A}{L} + \frac{q}{2k}(L - x) \right) x + T_A$

"""
#----------------------------------------------------------------------------------------
                               #Library Imports
#----------------------------------------------------------------------------------------

import FiniteVolumeMethod as fvm
import numpy as np
import matplotlib.pyplot as plt

#----------------------------------------------------------------------------------------
                         #Problem Variables Declaration
#----------------------------------------------------------------------------------------

longitud = 0.02                  # meters
TA = 100                         # °C 
TB = 200                         # °C 
k  = 0.5                         # W/m.K
q  = 1e+06                       # 1e+06 W/m^3 = 1000 kW/m^3 Fuente uniforme
N  = 6                           # Número de nodos

#----------------------------------------------------------------------------------------
                       #Mesh creation and mesh parameter definition
#----------------------------------------------------------------------------------------

malla = fvm.Mesh(nodes = N, length = longitud)
nx    = malla.nodes()                                # Número de nodos
nvx   = malla.volumes()                              # Número de volúmenes
delta = malla.delta()                                # Tamaño de los volúmenes

#----------------------------------------------------------------------------------------
                       #Problem variables printing
#----------------------------------------------------------------------------------------

fvm.printData(Longitud = longitud,
              Temperatura_A = TA,
              Temperatura_B = TB,
              Conductividad = k,
              Fuente = q,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)

#----------------------------------------------------------------------------------------
                #Diffusive coefficients allocation and calculation
#----------------------------------------------------------------------------------------

df1 = fvm.Diffusion1D(nvx, Gamma = k, dx = delta)
df1.alloc(nvx)                                      # Se aloja memoria para los coeficientes
df1.calcCoef()                                      # Se calculan los coeficientes
df1.source(q)                                       # Se agregua la fuente uniforme

#----------------------------------------------------------------------------------------
                    #Border conditions definition
#----------------------------------------------------------------------------------------

T = np.zeros(nvx)                       # El arreglo contiene ceros
T[0]  = TA                              # Condición de frontera izquierda
T[-1] = TB                              # Condición de frontera derecha
df1.bcDirichlet('LEFT_WALL', T[0])      # Se actualizan los coeficientes
df1.bcDirichlet('RIGHT_WALL', T[-1])    # de acuerdo a las cond. de frontera

#----------------------------------------------------------------------------------------
                        #Linear system solution matrix creation
#----------------------------------------------------------------------------------------

Su = df1.Su()                       # Vector del lado derecho
A = fvm.Matrix(malla.volumes())     # Matriz del sistema
A.build(df1)                        # Construcción de la matriz en la memoria

#----------------------------------------------------------------------------------------
                                 #Linear system solution
#----------------------------------------------------------------------------------------

T[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
print('Solución = ', T)

#----------------------------------------------------------------------------------------
                             #Graph x axis initialization
#----------------------------------------------------------------------------------------

x = malla.createMesh()

#----------------------------------------------------------------------------------------
                             #Analytycal Solution
#----------------------------------------------------------------------------------------

Ta = ((TB - TA) / longitud + q * (longitud - x) / (2 * k) ) * x + TA

#----------------------------------------------------------------------------------------
                              #Results Plotting
#----------------------------------------------------------------------------------------

x *= 100                                                    # Transformación a [cm]
plt.plot(x,Ta, '-', label = 'Sol. analítica')               # Sol. analítica
plt.plot(x,T,'o', label = 'Sol. FVM')
plt.title('Ecuación de Calor [example02]')
plt.ylim(75,300)
plt.xlabel('$x$ [cm]')
plt.ylabel('Temperatura [$^o$C]')
plt.grid()
plt.legend()
plt.savefig('example02.pdf')
plt.show()
