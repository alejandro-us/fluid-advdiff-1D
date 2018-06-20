#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat April 21 17:00:26 2018

@author: alex
"""
import numpy as np
from Coefficients import Coefficients

class Temporal1DBackEuler(Coefficients):
    """
    Class that defines coefficients originated from advection-convection equation using Finite Volume method with 
    implicit Euler approximation.
    
    Parameters:
        nvx: Number of volumes
        rho: Density
        dx: Size space step in a constant grid
        dt: Size time step in a constant grid
    """
    
    def __init__(self, nvx = None, rho = None, dx = None, dt = None, tsteps = None):
        super().__init__(nvx)
        self._nvx = nvx
        if rho:
            self._rho = rho
        else:
            self._rho = 1 
        if dx:
            self._dx = dx
        else:
            self._dx = 0.1
        if dt:
            self._dt = dt
        else:
            self._dt = 0.1 
        if tsteps:
            self._tsteps = tsteps
        else:
            self._tsteps = 20

    def __del__(self):
        del(self._nvx)
        del(self._rho)
        del(self._dx)
        del(self._dt)
        del(self._tsteps)

    def dx(self):
         return self._dx

    def dt(self):
        return self._dt
    
    def rho(self):
        return self.rho
    
    def calcCoef(self, priorSol = None):
        """
        Calculates the proper coefficients for advective-Convective movement in backward Euler approximation.
        
        Return:
            None
        """
        aP = self.aP()
        Su = self.Su()
        rho = self._rho
        dx = self._dx
        dt = self._dt

        for i in range(1,self._nvx-1):
           aP[i] += rho*(dx/dt)
           Su[i] += priorSol[i]*(dx/dt)

