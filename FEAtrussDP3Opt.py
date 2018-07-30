#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 12:05:58 2018
FEA class for trusses
@author: nikorose
"""

import numpy as np
import vtk
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from latex_envs.latex_envs import figcaption

class Vizualization():

    def __init__(self):
        pass
    def displayTruss(self, elemNodes,nodeCords, stress=None, Area=None, \
                     name="Quantity", displace=False, displacements=None, \
                     scale = 1.0):
        '''
        Function for calculating stiffness matrix, displacements, stresses of 2D trusses
        with the VTK libraries. Original code proposed by SukhbinderSingh.com    
        '''
        if displace == True:
            # If statement for displaying deformations only 
            # converting displacements into a 2D array
            #Adding the displacement to the coords
            if displacements.all():
                raise ValueError('No displacement found, did you forget stated them?')
            else:
                nodeCords =  np.reshape(displacements, (-1, 2)) * scale + nodeCords ## Note that the scale is put here 
            
        pts = vtk.vtkPoints()
     
        for x,y in nodeCords:
            pts.InsertNextPoint(x,y,0.0)
     
        lines=vtk.vtkCellArray()
        for ii,jj in elemNodes:
            lines.InsertNextCell(2)
            lines.InsertCellPoint(ii)
            lines.InsertCellPoint(jj)
   
        stdata = vtk.vtkDoubleArray()
        stdata.SetName(name)
        for val in stress:
            stdata.InsertNextValue(val)
      
        grid = vtk.vtkPolyData()
        grid.SetPoints(pts)
        grid.SetLines(lines)
     
        grid.GetCellData().SetScalars(stdata)
      
     
        mapper=vtk.vtkPolyDataMapper()
        mapper.SetInputData(grid)
        mapper.SetScalarRange(np.min(stress),np.max(stress))

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
     
        sbar = vtk.vtkScalarBarActor()
        sbar.SetLookupTable(mapper.GetLookupTable())
        sbar.SetTitle(name)
     
        ren = vtk.vtkRenderer()
     
        ren.AddActor2D(sbar)
        ren.AddActor(actor)
     
        renwin= vtk.vtkRenderWindow()
        renwin.AddRenderer(ren)
        renwin.SetSize(900,500)
      
        iren=vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renwin)
        iren.Initialize()
        renwin.Render()
        iren.Start()
        
    def visualizeAreas(self, elementNodes, nodes, AreaOpt, supports, force, scale, stress, linewidth, maxiter):
        '''
        A function to plot the results in terms of areas (width of the line) and in terms
        of Stresses (Blue is compression, red is tension).
        '''
        #Normalize stresses from 0 to 1
        Stress_norm = np.concatenate((stress-np.min(stress))/(np.max(stress)-np.min(stress)))
        colors = pl.cm.jet(Stress_norm)
        plt.figure()
        count = 0
        len_ele = elementNodes[2,0,1] - elementNodes[2,0,0]
        Areaviz =  (AreaOpt-np.min(AreaOpt))/(np.max(AreaOpt)-np.min(AreaOpt)) * 5.0 +0.1
        figcaption('Topology Opt. of a 2D truss Structure at {} iterations.'.format(maxiter), label="fig:{}".format(maxiter))
        for i in elementNodes:
            plt.plot(i[:,0], i[:,1], color= colors[count], linewidth = Areaviz[count])
            count += 1
#        plt.plot(nodes[force,0], nodes[force,1], marker="v", markersize=8, color="red")
        for e in supports:
            plt.plot(nodes[e,0], nodes[e,1], marker="*", markersize=8, color="blue")
        plt.axis('off')
        plt.annotate('force', xy=(nodes[force,0], nodes[force,1]), xytext=(nodes[force,0], nodes[force,1]+len_ele/2.0),
            arrowprops=dict(facecolor='red', shrink=0.05))
#        cbar = plt.colorbar(stress, ticks=[-1, 0, 1])
#        cbar.ax.set_yticklabels(['< -1', '0', '> 1'])# vertically oriented colorbar
        #plt.title('Topology Opt. of a 2D truss Structure at {} iterations.'.format(maxiter))
        plt.show()
        return 

class FEA():    
    '''
    Contains all the functions to perform a Finite element 
    Analysis of 2D trusses
    '''
    def structure2(self, Numx, Numy, Coord = 10.0):
        '''
        Build the element array, containing the coordinates in each element
        '''
        if Numx == 0 or Numy == 0:
            raise ValueError('The size of the matrix should be greater than zero')
        else:
            #Coordinates Matrix Building
            x = np.linspace(0, Numx*Coord, Numx+1)
            y = np.linspace(0, Numy*Coord, Numy+1)
            xnodes, ynodes = np.meshgrid(x, y, indexing= 'ij')
            nodes = np.flip(np.array(np.meshgrid(y, x, indexing= 'ij')).T.reshape(-1,2),1) # We'll always have two dims in this problem
            # Flip is for doing a mirror (I dont't it worked in that way) of the array
            coordinate_grid = np.array([xnodes, ynodes])
            Ny = 0
            GeneralElem = []
            while Ny <= Numy-1:
                Nx = 0
                GeneralElem += [[coordinate_grid[:,Nx,Ny], coordinate_grid[:,Nx,Ny+1]]] # elem 0 in the Y border
                while Nx <= Numx-1:
                    GeneralElem += [[coordinate_grid[:,Nx,Ny], coordinate_grid[:,Nx+1,Ny+1]]] #element 1
                    GeneralElem += [[coordinate_grid[:,Nx,Ny+1], coordinate_grid[:,Nx+1,Ny]]] # elem 2
                    GeneralElem += [[coordinate_grid[:,Nx,Ny+1], coordinate_grid[:,Nx+1,Ny+1]]] # elem 3
                    GeneralElem += [[coordinate_grid[:,Nx+1,Ny+1], coordinate_grid[:,Nx+1,Ny]]] # elem 4
                    if Ny == 0:
                        GeneralElem += [[coordinate_grid[:,Nx,Ny], coordinate_grid[:,Nx+1,Ny]]] # elem 5 in the borders
                    Nx += 1
                Ny += 1       
        return np.array(GeneralElem), nodes
    
    def indexation(self, elementNodes, nodes):
        '''Returns the element matrix but with the node indexation'''
        index = []
        for e in elementNodes:
            # Finding the common index in this search
            i11 = np.where(nodes[:,0] == e[0,0])
            i12 = np.where(nodes[:,1] == e[0,1])
            i21 = np.where(nodes[:,0] == e[1,0])
            i22 = np.where(nodes[:,1] == e[1,1])
            i1 = np.intersect1d(i11,i12)[0]
            i2 = np.intersect1d(i21,i22)[0]
            index.append([i1,i2])
        return np.array(index)
        
    def solution(self, GDof,prescribedDof,stiffness, force):
        '''function to find solution in terms of global displacements 
        by AJM Ferreira on his book FEA with Matlab
        '''
        activeDof=np.setdiff1d(np.arange(GDof),prescribedDof)
        U=np.linalg.solve(stiffness[np.ix_(activeDof,activeDof)],force[np.ix_(activeDof)]);
        displacements=np.zeros((GDof,1))
        displacements[np.ix_(activeDof)]= U
        return displacements
    
    def formStiffness2Dtruss(self,A, GDof, elementNodes,nodes, indices, E):
        '''
        The code was inspired by the book MATLAB 
        Codes for Finite Element Analysis by
        AJM Ferreira, Springer
        '''
        stiffness=np.zeros((GDof,GDof))
        count = 0
        for e in elementNodes:
            elemDof=np.array([indices[count,0]*2, indices[count,0]*2+1, indices[count,1]*2, indices[count,1]*2+1])
            count += 1
            # Take into account the order of the next substraction
            xa = e[0,0]-e[1,0]
            ya = e[0,1]-e[1,1]
            len_elem=np.sqrt(xa**2+ya**2)
            c=xa/len_elem
            s=ya/len_elem
            try:
                EA = E*A[0,count-1]
            except IndexError:
                EA = E*A[count-1] ## The optimization algorithm change np array to a list
            k1=(EA/len_elem)* np.array([[c*c,c*s,-c*c, -c*s],
                                        [c*s,s*s,-c*s ,-s*s],
                                        [-c*c,-c*s,c*c,c*s],
                                        [-c*s,-s*s,c*s,s*s]])
            stiffness[np.ix_(elemDof,elemDof)] +=k1
        return stiffness
    
    def outputDisplacementsReactions(self, displacements,stiffness,GDof,prescribedDof):
        '''
        Calculates the force reactions once the displacements are obtained
        '''
        F = np.dot(stiffness,displacements)
        reactions = F[prescribedDof.T][0]
#        print(np.hstack([prescribedDof, reactions[:]]))
        return reactions
    
    def stresses2Dtruss(self, elementNodes,indices, displacements,E):
        '''
        Calculates the stresses in each element
        '''
        sigma = np.zeros((elementNodes.shape[0],1))
        count = 0
        for e in elementNodes:
            elemDof=np.array([indices[count,0]*2, indices[count,0]*2+1, indices[count,1]*2, indices[count,1]*2+1])
#            xa=np.abs(e[0,0]-e[1,0])
#            ya=np.abs(e[0,1]-e[1,1])
            xa = e[1,0]-e[0,0]
            ya = e[1,1]-e[0,1]
            len_elem=np.sqrt(xa**2+ya**2)
            c=xa/len_elem
            s=ya/len_elem
            sigma[count] = (E/len_elem) * np.dot(np.array([-c,-s,c,s]),displacements[np.ix_(elemDof)])
            count += 1
        return sigma
    

