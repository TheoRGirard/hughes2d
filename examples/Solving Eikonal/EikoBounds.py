from hughes2d import *
import random as alea
import numpy as np
import plotly.graph_objects as go

"""MyDomain = Mesh2D.NonConvexDomain()
MyDomain = Mesh2D.NonConvexDomain([[0,0],[10,0],[10,3],[12,3],[12,4],[10,4],[10,7],[0,7]])
MyDomain.addExit([[12,3],[12,4]])

MyMesh = Mesh2D.Mesh()
MyMesh.generateMeshFromDomain(MyDomain, 0.015)
MyMesh.show()
"""
#MyMesh.saveToJson("test1")

MyMesh = Mesh2D.Mesh()

MyMesh.loadFromJson("test1_mesh.json")

#MyMesh.show()

MyMap = Mesh2D.CellValueMap(MyMesh)
#MyMap.setConstantCircle([7,4],2.5,0.6)
#MyMap.show()

ConstantVectorField = []

for P in MyMesh.barycenters:
    if(3<=P[1]<=4):
        ConstantVectorField.append([1,0])
    elif(P[1]>4):
        V = [10-P[0],4-P[1]]
        ConstantVectorField.append([V[0]/np.sqrt(V[0]**2 + V[1]**2), V[1]/np.sqrt(V[0]**2 + V[1]**2)])
    elif(P[1]<3):
        V = [10-P[0],3-P[1]]
        ConstantVectorField.append([V[0]/np.sqrt(V[0]**2 + V[1]**2), V[1]/np.sqrt(V[0]**2 + V[1]**2)])

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
