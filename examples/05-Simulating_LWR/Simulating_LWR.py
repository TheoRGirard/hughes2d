import random as alea

import matplotlib.pyplot as plt
import numpy as np
import plotly

from hughes2d import *

#Parameters of the simulation
dt = 0.025
dx = 0.05
numStep = 200

#Construction of the mesh
MyDomain = Mesh2D.NonConvexDomain([[0,0],[10,0],[10,5],[0,5]])
MyDomain.addExits([[[10,0],[10,5]]])
MyMesh = Mesh2D.Mesh()
MyMesh.generateMeshFromDomain(MyDomain,dx)

#Construction of an initial datum
InitialDatum = Mesh2D.CellValueMap(MyMesh)
InitialDatum.generateRandom()

#Plot the InitialDatum
InitialDatum.show(preference="matplotlib")

#Construction of a vector field
VectorField = [[1,0] for _ in MyMesh.triangles]

#Construction of the LWR solver
opt = dict( method = "midVector",
            convexFlux = True )
MySolver = LWR2D.LWRSolver(MyMesh, dt, previousDensity=InitialDatum, directionMap =VectorField, opt = opt)

#Compute the number of steps numStep
for _ in range(numStep):
    MySolver.computeNextStep()
    MySolver.update(VectorField)

#Plot the final state of the density
MySolver.showDensity(t=1, preference="matplotlib")
