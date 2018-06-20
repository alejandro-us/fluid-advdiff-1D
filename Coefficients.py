#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 15:11:05 2018

@author: luiggi
modified by alex
"""

import numpy as np

class Coefficients():
    """
    This class defines the main arrays needed for the Finite Volume Method.
    The arraysare defined as class variables to be shared with all instances 
    of such class.
    """    
    _aP = None
    _aE = None
    _aW = None
    _aEE = None
    _aWW = None
    _Su = None
    _Sp = None
    _nvx = None
    _delta = None

    def __init__(self, nvx = None, delta = None):
        """
        Class Constructor
        
         Param:
            nvx: number of volumes
            delta: constant grid distance step lenght
        Return:
            None
        """
        Coefficients._nvx = nvx
        Coefficients._delta = delta


    @staticmethod
    def alloc(n):
        """
        allocate space in memory for coeficients.
        
        Param:
            n: number of volumes
        Return:
            None
            
        """
        if Coefficients._nvx:
            nvx = Coefficients._nvx
        else:
            nvx = n
        Coefficients._aP = np.zeros(nvx)
        Coefficients._aE = np.zeros(nvx)
        Coefficients._aW = np.zeros(nvx)
        Coefficients._Su = np.zeros(nvx)
        Coefficients._Sp = np.zeros(nvx)
        Coefficients._aEE = np.zeros(nvx)
        Coefficients._aWW = np.zeros(nvx)
    
    def setVolumes(self, nvx):
        Coefficients._nvx = nvx
        
    def setDelta(self, delta):
        Coefficients._delta = delta
        
    def aP(self):
        return Coefficients._aP

    def aE(self):
        return Coefficients._aE
    
    def aW(self):
        return Coefficients._aW
    
    def Su(self):
        return Coefficients._Su
    
    def Sp(self):
        return Coefficients._Sp
    
    def aEE(self):
        return Coefficients._aEE
    
    def aWW(self):
        return Coefficients._aWW
    
    def cleanCoefficients(self):
        """
        Assigns zero to all coefficients.
        
        Param:
            None
        Return:
            None
        """
        Coefficients._aP[:] = 0.0
        Coefficients._aE[:] = 0.0
        Coefficients._aW[:] = 0.0
        Coefficients._Su[:] = 0.0

    @staticmethod
    def bcDirichlet(wall, phi):
        """
        Assign values to the boundary of the grid modifying coefficients.
        
        Param:
            wall: literal used to decide which boundary will be modified
            phi: value of dirichlet boundarie value
        Return:
            None
        """
        aP = Coefficients._aP
        aE = Coefficients._aE
        aW = Coefficients._aW
        Su = Coefficients._Su

        if wall == 'LEFT_WALL':
            aP[1] += aW[1]
            Su[1] += 2 * aW[1] * phi
        elif wall == 'RIGHT_WALL':
            aP[-2] += aE[-2]
            Su[-2] += 2 * aE[-2] * phi    
            
            
    @staticmethod
    def bcNeumann(wall, phi):
        """
        Assign values to the boundary of the grid modifying coefficients.
        
        Param:
            wall: literal used to decide which boundary will be modified
            phi: value of Neumann boundary value
        Return:
            None
        """
        aP = Coefficients._aP
        aE = Coefficients._aE
        aW = Coefficients._aW
        Su = Coefficients._Su
        delta = Coefficients._delta

        if wall == 'LEFT_WALL':
            aP[1] -= aW[1]
            Su[1] += delta * aW[1] * phi
        elif wall == 'RIGHT_WALL':
            aP[-2] -= aE[-2]
            Su[-2] += delta * aE[-2] * phi    
    
    @staticmethod   
    def source(q = None, s = None):
        """
        Defines the source values directly across the grid domain.
        
        Param:
            q: fix source value
            s: linear source component (source with property value dependance)
        Return:
            None
        """
        Sp = Coefficients._Sp
        Su = Coefficients._Su
        dx = Coefficients._delta
        Su += q * dx
        if s:
            Sp += s 
        

if __name__ == '__main__':
    
    coef1 = Coefficients(6, 0.25)
    coef1.alloc(6)
    coef1.source(100)
    
    print('-' * 20)  
    print(coef1.aP(), coef1.aE(), coef1.aW(), coef1.Su(), sep = '\n')
    print('-' * 20)  

    ap = coef1.aP()
    ap[2] = 25
    print(ap, coef1.aP(),sep='\n')
    print('-' * 20)  

    ae = coef1.aE()
    aw = coef1.aW()
    su = coef1.Su()
    ae.fill(5)
    aw.fill(5)
    ap.fill(10)
    coef1.bcDirichlet('LEFT_WALL', 2)
    coef1.bcDirichlet('RIGHT_WALL', 1)
    print(coef1.aP(), coef1.aE(), coef1.aW(), coef1.Su(), sep = '\n')
    print('-' * 20)  

