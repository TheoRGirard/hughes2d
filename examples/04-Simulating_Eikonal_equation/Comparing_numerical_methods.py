import random as alea

import numpy as np
import plotly.graph_objects as go

from hughes2d import *

#Loading the mesh ---------------------------
MyMesh = Mesh2D.Mesh()
MyMesh.loadFromJson("examples/04-Simulating_Eikonal_equation/data/rectangular_mesh.json")
MyEmptyMap = Mesh2D.CellValueMap(MyMesh)

#Constructing the explicit solution
ExplicitSolution = []
for P in MyMesh.vertices:
    if(3<=P[1]<=4):
        ExplicitSolution.append(12-P[0])
    elif(P[1]>4):
        V = [10-P[0],4-P[1]]
        ExplicitSolution.append(np.sqrt(V[0]**2 + V[1]**2)+2)
    elif(P[1]<3):
        V = [10-P[0],3-P[1]]
        ExplicitSolution.append(np.sqrt(V[0]**2 + V[1]**2)+2)

MySol = Mesh2D.VertexValueMap(MyMesh)
MySol.values = ExplicitSolution

Opt=dict(method = "FMT", constrained = False, NarrowBandDepth = 2)

FMTuSolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt)
FMTuSolv.computeField()

Opt3=dict(method = "FMT", constrained = True, NarrowBandDepth = 2)

FMTCSolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt3)
FMTCSolv.computeField()

Opt2=dict(method = "FME", constrained = True, NarrowBandDepth = 2)

FMEsolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt2)
FMEsolv.computeField()


fig = go.Figure()

FMTuSolv.fieldValues.add3Dplot(fig,color=[255,127,0,0.5])
FMTCSolv.fieldValues.add3Dplot(fig,color=[0,255,0,0.5])
FMEsolv.fieldValues.add3Dplot(fig,color=[0,0,255,0.5])

MySol.add3Dplot(fig,color=[255,0,0,0.5])

fig.show()
