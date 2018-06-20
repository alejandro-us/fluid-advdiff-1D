#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 13:21:50 2018

@author: luiggi
modified by alex
"""

import numpy as np

class Mesh():
    """
        Class thet contains spatial information about the greed including number of nodes and volumes.
           
        Atributes:
            
            Nodes: Number of nodes in de 1 dimensional mesh.
            Volumes: Number of volumes in the 1 dimensional mesh.
            Lenght: Lenght of the mesh.
            Delta: Lenght of the disntance within volumes of a equally distributed mesh.
            timeInte: Time interval used in temporal solution for advective-diffusive problem.
            tsteps: Number of time steps (iterations) used to solve temporal advective-diffusive problem.
    
    """
    
    def __init__(self, nodes = None, volumes = None, length = None, timeInt = None, tsteps = None ):
        """
            Class Constructor

        """
        self.__nodes = nodes
        self.__volumes = volumes
        self.__length = length     
        self.__delta = 1
        self.adjustNodesVolumes(nodes, volumes)
        self.calcDelta()
        self.__timeInt = timeInt
        self.__tsteps = tsteps
    
    def __del__(self):
        
        """
            Class Destructor
         
        """
        del(self.__nodes)
        del(self.__volumes)
        del(self.__length)
        del(self.__delta)
        del(self.__timeInt)
        del(self.__tsteps)
        
    def adjustNodesVolumes(self,nodes,volumes):
        if nodes:
            self.__volumes = self.__nodes + 1
        if volumes:
            self.__nodes = self.__volumes - 1        
        
    def nodes(self):
        return self.__nodes
    
    def setNodes(self, nodes):
        self.__nodes = nodes
        self.adjustNodesVolumes(nodes = nodes, volumes = None)
        
    def volumes(self):
        return self.__volumes

    def setVolumes(self, volumes):
        self.__volumes = volumes
        self.adjustNodesVolumes(nodes = None, volumes = volumes)
        
    def length(self):
        return self.__length
        
    def calcDelta(self):
        if self.__length:
            self.__delta = self.__length / (self.__nodes - 1)
        
    def delta(self):
        return self.__delta
    
    def timeInt(self):
        return self.__timeInt
    
    def tsteps(self):
        return self.__tsteps
    
    def createMesh(self):
        '''
        Creates the data array with euclidian space inforation about the grid.
        '''
        first_volume = self.__delta / 2
        final_volume = self.__length - first_volume
        self.__x = np.zeros(self.__volumes)
        self.__x[1:-1] = np.linspace(first_volume,final_volume,self.__volumes-2)
        self.__x[-1] = self.__length
        return self.__x
    
    def createTimeMesh(self):
        '''
        Creates the data spline corresponding to the time interval and time steps of the temporal solution of the problem.
        '''
        self.__t = np.linspace(0,self.__timeInt, self.__tsteps + 1)
        return self.__t
        
if __name__ == '__main__':

    m1 = Mesh()
    print(m1.nodes(), m1.volumes())
    print('_' * 20)   
    
    m1 = Mesh(nodes = 5)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)
   
    m1 = Mesh(volumes = 5)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)
    
    m1 = Mesh(5,5)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)
    
    m1.setNodes(8)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)

    m1.setVolumes(8)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)
    
    m1 = Mesh(nodes =  5, length = 33)
    print(m1.nodes(), m1.volumes(), m1.length())
    print('_' * 20)
    
    m1 = Mesh(volumes =  5, length = 33)
    print(m1.nodes(), m1.volumes(), m1.length())
    print('_' * 20)
    
    m1 = Mesh(nodes = 5, length = 1)
    print(m1.nodes(), m1.volumes(), m1.length(), m1.delta())
    print('_' * 20)    
    
    m1 = Mesh(volumes = 10, length = 1)
    print(m1.nodes(), m1.volumes(), m1.length(), m1.delta())
    print('_' * 20) 

    m1 = Mesh(volumes = 6, length = 1)
    print(m1.nodes(), m1.volumes(), m1.length(), m1.delta())
    m1.createMesh()
    print('_' * 20) 
    
    m1 = Mesh(volumes = 6, length = 1, timeInt=1 , tsteps=500)
    print(m1.nodes(), m1.volumes(), m1.length(), m1.delta())
    print(m1.timeInt(), m1.tsteps())
    print(m1.createTimeMesh())
    print('_' * 20) 