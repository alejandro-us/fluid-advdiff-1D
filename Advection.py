#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 18:46:43 2018

@author: luiggi
modified by alex
"""

import numpy as np
from Coefficients import Coefficients

#                           nx = 5, nvx = 6
#                    0     1     2     3     4
#                    x--|--x--|--x--|--x--|--x     <-- Velocities dephased - 1/2 (Uw = u[i-1], Ue = u[i])
#                 0     1     2     3     4     5  <-- Volumes 
#                 o--|--x--|--x--|--x--|--x--|--o
#                       0     1     2     3        <-- Unknowns  

class Advection1DUpwind(Coefficients):
    """
    Class that defines coefficients originated from advection equation using Finite Volume method with 
    upwind approximation.
    
    Parameters:
        nvx: Number of volumes\n
        rho: Density\n
        dx: Size space step in a constant grid\n
        u: velocity
    """
    
    def __init__(self, nvx = None, rho = None, dx = None, u = None):
        super().__init__(nvx)
        self._nvx = nvx
        if rho:
            self._rho = rho
        else:
            self._rho = 1
        self._dx = dx
        if not (u is None):
            self._u = u
        else:
            self._u = np.zeros(nvx-1)

    def __del__(self):
        del(self._nvx)
        del(self._rho)
        del(self._dx)
        del(self._u)

    def setVel(self, u):
        self._u = u

    def u(self):
        return self._u
    
    def calcCoef(self):
        """
        Calculates the proper coefficients for advective movement in upwind approximation.
        
        Return:
            None
        """
        aE = self.aE()
        aW = self.aW()
        aP = self.aP()
        u = self._u
        rho = self._rho

        for i in range(1,self._nvx-1):
           aE[i] += max(-u[i]*rho,0)
           aW[i] += max(u[i-1]*rho,0)
           aP[i] += max(-u[i]*rho,0) + max(u[i-1]*rho,0) + rho * (u[i] - u[i-1])
  
#-------------------------------------------------------------------------------------------------------------------------          

class Advection1DCDS(Coefficients):
    """
    Class that defines coefficients originated from advection equation using Finite Volume method with 
    Central Difference Scheme (CDS) approximation.
    
    Parameters:
        nvx: Number of volumes
        rho: Density
        dx: Size space step in a constant grid
        u: velocity
    """
    
    def __init__(self, nvx = None, rho = None, dx = None, u = None):
        super().__init__(nvx)
        self._nvx = nvx
        if rho:
            self._rho = rho
        else:
            self._rho = 1
        self._dx = dx
        if not (u is None):
            self._u = u
        else:
            self._u = np.zeros(nvx-1)

    def __del__(self):
        del(self._nvx)
        del(self._rho)
        del(self._dx)
        del(self._u)

    def setVel(self, u):
        self._u = u

    def u(self):
        return self._u
    
    def calcCoef(self):
        """
        Calculates the proper coefficients for advective movement in CDS approximation.
        
        Return:
            None
        """
        aE = self.aE()
        aW = self.aW()
        aP = self.aP()
        u = self._u
        rho = self._rho

        for i in range(1,self._nvx-1):
            aE[i] += (-u[i]*rho)/2
            aW[i] += (u[i-1]*rho)/2
            aP[i] += (-u[i]*rho)/2 + (u[i-1]*rho)/2 + rho * (u[i] - u[i-1])

#----------------------------------------------------------------------------------------------------------
            
class Advection1DQuick(Coefficients):
    """
    Class that defines coefficients originated from advection equation using Finite Volume method with 
    quadratic upstream interpolation for convective kinetics (QUICK) approximation.
    
    Parameters:
        nvx: Number of volumes
        rho: Density
        dx: Size space step in a constant grid
        u: velocity
    """  
#                
    def __init__(self, nvx = None, rho = None, dx = None, u = None):
        super().__init__(nvx)
        self._nvx = nvx
        if rho:
            self._rho = rho
        else:
            self._rho = 1
        self._dx = dx
        
        if not (u is None):
            self._u = u
        else:
            self._u = np.zeros(nvx-1)

    def __del__(self):
        del(self._nvx)
        del(self._rho)
        del(self._dx)
        del(self._u)

    def setVel(self, u):
        self._u = u

    def u(self):
        return self._u
    
    def calcCoef(self):
        """
        Calculates the proper coefficients for advective movement in Quadratic upwind differencing
        scheme (QUICK)  approximation.
        
        Return:
            None
        """
        aE = self.aE()
        aW = self.aW()
        aP = self.aP()
        aEE = self.aEE()
        aWW = self.aWW()
        u = self._u
        rho = self._rho

        for i in range(3,self._nvx-2):
            
            ce = u[i]*rho
            cw = u[i-1]*rho
            
            if ce > 0:
                alfa_e = 1
            else:
                alfa_e = 0
            if cw > 0:
                alfa_w = 1
            else:
                alfa_w = 0
            
            Cw = 6/8*(alfa_w * cw) + (1/8)* alfa_e *ce + (3/8)*(1-alfa_w)*cw
            Cww = -(1/8)*alfa_w*cw
            Ce = -3/8*(alfa_e * ce) - (6/8)* (1-alfa_e)*ce - (1/8)*(1-alfa_w)*cw
            Cee = (1/8)*(1-alfa_e)*ce
            
            aE[i] += Ce
            aW[i] += Cw
            aEE[i] += Cee
            aWW[i] += Cww
            aP[i] += Ce + Cw + Cww + Cee + (ce - cw)

     
    def bcDirichlet(self, wall, phi):
        """
        Assign values to the boundary of the grid modifying coefficients using QUICK approximation scheme.
        (Overrride)
        Param:
            wall: literal used to decide which boundary will be modified
            phi: value of dirichlet boundarie value
        Return:
            None
        """
        aP = self.aP()
        aE = self.aE()
        aW = self.aW()
        Su = self.Su()
        aEE = self.aEE()
        aWW = self.aWW()
        u = self._u
        rho = self._rho
        Da = aW[0]
        Db = aW[-1]
        

        if wall == 'LEFT_WALL':
            #First Node-----------------------------------------------
            aWW[1] = 0
            aW[1] = 0
            aE[1] += (1/3)*Da - (3/8)*rho*u[1]
            aP[1] = aW[1] + aE[1] + rho*(u[1]-u[0]) + (8/3)*Da +(2/8)*rho*u[1] + rho*u[0]
            Su[1] += ((8/3)*Da + (2/8)*rho*u[1] + rho*u[0])*phi
            #Second Node-----------------------------------------------
            aWW[2] = 0
            aW[2] += (7/8)*rho*u[1] + (1/8)*rho*u[2]
            aE[2] += -(3/8)*rho*u[2]
            aP[2] = aW[2] + aE[2] + rho*(u[2]-u[1]) - (1/4)*rho*u[1]
            Su[2] += -(1/4)*rho*u[1]*phi
            
        elif wall == 'RIGHT_WALL':
            #Last Node----------------------------------------------------
            aWW[-2] = -(1/8)*rho*u[-2]
            aW[-2] += (1/3)*Db + (6/8)*rho*u[-2]
            aE[-2] = 0
            aP[-2] = aWW[-2] + aW[-2] + aE[-2] + rho*(u[-2]-u[-1]) + (8/3)*Db - rho*u[-1]
            Su[-2] += ((8/3)*Db - rho*u[-1])*phi
#--------------------------------------------------------------------------------------------------------

class Advection1DUpwind2(Coefficients):
    """
    Class that defines coefficients originated from advection equation using Finite Volume method with 
    Second order Upwind approximation.
    
    Parameters:
        nvx: Number of volumes
        rho: Density
        dx: Size space step in a constant grid
        u: velocity
    """
    
    def __init__(self, nvx = None, rho = None, dx = None, u = None):
        super().__init__(nvx)
        self._nvx = nvx
        if rho:
            self._rho = rho
        else:
            self._rho = 1
        self._dx = dx
        if not (u is None):
            self._u = u
        else:
            self._u = np.zeros(nvx-1)

    def __del__(self):
        del(self._nvx)
        del(self._rho)
        del(self._dx)
        del(self._u)

    def setVel(self, u):
        self._u = u

    def u(self):
        return self._u
    
    def calcCoef(self):
        """
        Calculates the proper coefficients for advective movement in CDS approximation.
        
        Return:
            None
        """
        aE = self.aE()
        aW = self.aW()
        aP = self.aP()
        aEE = self.aEE()
        aWW = self.aWW()
        u = self._u
        rho = self._rho

        for i in range(1,self._nvx-1):
            
            ce = u[i]*rho
            cw = u[i-1]*rho
            
            aE[i] += (-3/2)*max((-ce,0)) + (-1/2)*max(-cw,0)
            aW[i] += (1/2)*max((ce,0)) + (3/2)*max((cw,0))
            aEE[i] += (1/2)*max((-ce,0))
            aWW[i] += (-1/2)*max((cw,0))
            aP[i] += (3/2)*(max((ce,0)) + (-3/2)*max((-cw,0)))
            
    def bcDirichlet(self, wall, phi):
        """
        Assign values to the boundary of the grid modifying coefficients using Second order Upwind approximation scheme.
        (Overrride)
        Param:
            wall: literal used to decide which boundary will be modified
            phi: value of dirichlet boundarie value
        Return:
            None
        """
        aP = self.aP()
        aE = self.aE()
        aW = self.aW()
        Su = self.Su()
        aEE = self.aEE()
        aWW = self.aWW()
        u = self._u
        rho = self._rho  

        if wall == 'LEFT_WALL':
            #First Node-----------------------------------------------
               aP[1] += aW[1] + 3*aWW[1]
               Su[1] += (2*aW[1] + 4*aWW[1])*phi
               
            #Second Node-----------------------------------------------
               aW[2] += -aWW[2]
               Su[2] += 2*aWW[2]*phi
            
        elif wall == 'RIGHT_WALL':
             #Last Node-----------------------------------------------
               aP[-2] += aE[-2] + 3*aEE[-2]
               Su[-2] += (2*aE[-2] + 4*aEE[-2])*phi
            
            #Second lat Node-----------------------------------------------
               aE[-3] += -aEE[-3]
               Su[-3] += 2*aEE[-3]*phi
            
#--------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    
    nx = 5
#    u = np.sin(np.linspace(0,1,nx))
    u = np.ones(nx)
    print('-' * 20)  
    print(u)
    print('-' * 20)  

    af1 = Advection1DUpwind(6, 1, 1)
    af1.alloc(6)
    af1.setVel(u)
    print(af1.u())
    print('-' * 20)  

    af1.calcCoef()
    print(af1.aP(), af1.aE(), af1.aW(), af1.Su(), sep = '\n')
    print('-' * 20)  

    af1.bcDirichlet('LEFT_WALL', 2)
    af1.bcDirichlet('RIGHT_WALL', 1)
    print(af1.aP(), af1.aE(), af1.aW(), af1.Su(), sep = '\n')
    print('-' * 20)  



