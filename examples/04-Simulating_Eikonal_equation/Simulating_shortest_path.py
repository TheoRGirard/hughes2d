from hughes2d import *
import random as alea
import numpy as np

MyMesh = Mesh2D.Mesh()

#Construction of the mesh
#MyDomain = Mesh2D.NonConvexDomain([[0,0],[3,0],[3,3],[0,3]])
#MyDomain.addExits([[[1,0],[2,0]],[[3,0.5],[3,1.5]],[[0.5,3],[2,3]]])
#MyMesh.generateMeshFromDomain(MyDomain, 0.01)
#MyMesh.saveToJson("examples/04-Simulating_Eikonal_equation/data/square")

#Loading the mesh
MyMesh.loadFromJson("examples/04-Simulating_Eikonal_equation/data/square_mesh.json")

MyMap = Mesh2D.CellValueMap(MyMesh)

Opt=dict(method = "FMT", constrained = True, NarrowBandDepth = 2)

EikoSolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt)
EikoSolv.computeField()

#Show as a colored scatter plot
EikoSolv.fieldValues.show(grid=True, colorscale_name = 'magenta')

#Show as a 3D plot
EikoSolv.fieldValues.show3D()

#Displays the vector field associated with the computed approximation
EikoSolv.fieldValues.showVectorField()
