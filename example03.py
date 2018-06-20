#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 15:28:39 2018

@author: alex
Example 4.3 from Malalasekera Book
----------------------------------
Convection gives rise to a temperature-dependent heat loss or sink term in the 
governing equation. Shown in Figure 4.9 is a cylindrical fin with uniform cross-
sectional area A. The base is at a temperature of 100°C (T B ) and the end is
insulated. The fin is exposed to an ambient temperature of 20°C. One-dimensional 
heat transfer in this situation is governed by:

\frac{d}{dx}\left ( kA\frac{dT}{dx} \right )- hP\left ( T-T_{\infty } \right )= 0


where h is the convective heat transfer coefficient, P the perimeter, k the
thermal conductivity of the material and T ∞ the ambient temperature.
Calculate the temperature distribution along the fin and compare the results
with the analytical solution given by:

\frac{T-T_{\infty }}{T_{B}-T_{\infty }} = \frac{cosh\left [ n \left ( L-x \right )\right ]}{cosh\left ( nL \right )}


where n 2 = hP/(kA), L is the length of the fin and x the distance along the
fin. Data: L = 1 m, hP/(kA) = 25/m 2 (note that kA is constant).

___________________________________________________________

 |<-------------- 1.0 m ------------->|
 B                                    
 |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
 |                                    |  dT/Dx = 0
 |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

T_B = 100                            T_amb = 20
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

longitud = 1             # meters
Tamb = 20                # °C 
TB = 100                 # °C 
n2  = 25                 # hP/(kA) = 25 m*2
N  = 6                   # Número de nodos
k = 1                    #la gamma es constante en este ejemplo ya que es factorizada para calcular n2
q = Tamb*n2              #la fuente en el centro de la malla
DtDx = 0                 #la razon de cambio en la frontera derecha

#----------------------------------------------------------------------------------------
                       #Mesh creation and mesh parameter definition
#----------------------------------------------------------------------------------------

malla = fvm.Mesh(nodes = N, length = longitud)
nx    = malla.nodes()                               # Número de nodos
nvx   = malla.volumes()                             # Número de volúmenes
delta = malla.delta()                               # Tamaño de los volúmenes

#----------------------------------------------------------------------------------------
                       #Problem variables printing
#----------------------------------------------------------------------------------------

fvm.printData(Longitud = longitud,
              Temperatura_ambiente = Tamb,
              Temperatura_B = TB,
              Condicion_Neumann = DtDx,
              Conductividad = k,
              Fuente = q,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)

#----------------------------------------------------------------------------------------
                #Diffusive coefficients allocation and calculation
#----------------------------------------------------------------------------------------

df1 = fvm.Diffusion1D(nvx, Gamma = k, dx = delta)
df1.alloc(nvx)                                          # Se aloja memoria para los coeficientes
Coefficients.source(q,-n2*delta)                        # Se agrega la fuente uniforme
df1.calcCoef()                                          # Se calculan los coeficientes


#----------------------------------------------------------------------------------------
                    #Border conditions definition
#----------------------------------------------------------------------------------------

T = np.zeros(nvx)                           # El arreglo contiene ceros
T[0]  = TB                                  # Condición de frontera izquierda
T[-1] = Tamb                                # Aproximacion de la condicion de frontera derecha usando la condicion tipo Neumman

df1.bcDirichlet('LEFT_WALL', T[0])          # Se actualizan los coeficientes
df1.bcNeumann('RIGHT_WALL', DtDx)           # de acuerdo a las cond. de frontera

#----------------------------------------------------------------------------------------
                        #Linear system solution matrix creation
#----------------------------------------------------------------------------------------

Su = df1.Su()                           # Vector del lado derecho
A = fvm.Matrix(malla.volumes())         # Matriz del sistema
A.build(df1)                            # Construcción de la matriz en la memoria

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

Ta = (np.cosh(((n2)**.5)*(longitud - x))/np.cosh((n2**.5)*longitud))*(TB-Tamb)+Tamb

#----------------------------------------------------------------------------------------
                              #Results Plotting
#----------------------------------------------------------------------------------------

x *= 100                                                        # Transformación a [cm]
plt.plot(x,Ta, '-', label = 'Sol. analítica')                   # Sol. analítica
plt.plot(x,T,'o', label = 'Sol. FVM')
plt.title('Ecuación de Calor [example03]')
plt.xlabel('$x$ [cm]')
plt.ylabel('Temperatura [$^o$C]')
plt.grid()
plt.legend()
plt.savefig('example03.pdf')
plt.show()


