#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 18:15:28 2018

@author: luiggi
modified by alex
"""
import numpy as np

class Matrix():
    """
    Creates and object in the form of a two dimensional array with the values of FVM coefficients.
    
    Atributes:
        N: size of the matrix
        A: Two dimensional array
    """
    
    def __init__(self, nvx = None):
        """
        Class Cosntructor
        """
        self.__N = nvx - 2 
        self.__A = np.eye(self.__N)

    def __del__(self):
        """
        Class Destructor
        """
        del(self.__N)
        del(self.__A)
        
    def mat(self):
        """
        Returns the array with matrix values.
        """
        return self.__A
    
    def build(self, coefficients = None):
        """
            Fill the matrix array with the corresponding coeficcients in order to solve linear equation system.
            
            Param:
                coefficientes: np.array
                
            Returns:
                None
         """
#                 nx = 5, nvx = 6
#                 0     1     2     3     4     5  <-- Volumes 
#                 o--|--x--|--x--|--x--|--x--|--o
#                       0     1     2     3        <-- Unknowns    
#                

        aP = coefficients.aP()
        aE = coefficients.aE()
        aW = coefficients.aW()
        aWW = coefficients.aWW()
        aEE = coefficients.aEE()
        A = self.__A
        A[0][0] = aP[1]
        A[0][1] = -aE[1]
        
        for i in range(1,self.__N-1): # range(1,N-3)  <-- (1,2)
            A[i][i] = aP[i+1]
            if i-2 >= 0:
                A[i][i-2] = -aWW[i+1]
            A[i][i+1] = -aE[i+1]
            A[i][i-1] = -aW[i+1]
            if i+2 <= self.__N-1:
                A[i][i+2] = -aEE[i+1]
                
        A[-1][-1] = aP[-2]
        A[-1][-2] = -aW[-2]

if __name__ == '__main__':

    a = Matrix(6)
    print('-' * 20)  
    print(a.mat())
    print('-' * 20)  
    
    from Diffusion import Diffusion1D        
    df1 = Diffusion1D(6, 1, 0.25)
    df1.alloc(6)
    df1.calcCoef()
    df1.source(100)
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  

    df1.bcDirichlet('LEFT_WALL', 2)
    df1.bcDirichlet('RIGHT_WALL', 1)
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  
    
    a.build(df1)
    print(a.mat())
    print('-' * 20)  
