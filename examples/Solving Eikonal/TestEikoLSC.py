from hughes2d import *
import random as alea
import numpy as np
import plotly.graph_objects as go
"""
MyDomain = Mesh2D.NonConvexDomain()
MyDomain = Mesh2D.NonConvexDomain([[-1,0],[-1,2],[1,2],[1,0]])
MyDomain.addExit([[-1,0],[1,0]])

MyMesh = Mesh2D.Mesh()
MyMesh.generateMeshFromDomain(MyDomain, 0.002)
MyMesh.show()

MyMesh.saveToJson("test2")
"""
MyMesh = Mesh2D.Mesh()

MyMesh.loadFromJson("test2_mesh.json")

#MyMesh.show()

MyMap = Mesh2D.CellValueMap(MyMesh)

for index, triangle in enumerate(MyMesh.triangles):
    MyMap[index] = 1
    for pair in [[0,1],[1,2],[2,0]]:
        if(MyMesh.vertices[triangle[pair[0]]][0])*(MyMesh.vertices[triangle[pair[1]]][0]) <= 0:
            MyMap.values[index] = 0
#MyMap.show()


ConstantVectorField = []

for P in MyMesh.barycenters:
    if(3*abs(P[0])/np.sqrt(3)>=P[1]):
        ConstantVectorField.append([0,-1])
    elif(abs(P[0]) < 0.001):
        ConstantVectorField.append([0,-1])
    else:
        V = [-P[0],abs(P[0]/3)]
        ConstantVectorField.append([V[0]/np.sqrt(V[0]**2 + V[1]**2), V[1]/np.sqrt(V[0]**2 + V[1]**2)])

ExplicitSolution = []

for P in MyMesh.vertices:
    if(3*abs(P[0])/np.sqrt(3)>=P[1]):
        ExplicitSolution.append(2*P[1])
    elif(abs(P[0]) < 0.0001):
        ExplicitSolution.append(P[1])
    else:
        ExplicitSolution.append(P[1] + 3*abs(P[0])/np.sqrt(3))


MySol = Mesh2D.VertexValueMap(MyMesh)
MySol.values = ExplicitSolution

Opt=dict(method = "FMT", constrained = True, NarrowBandDepth = 1)

EikoSolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt)
#EikoSolv.showNarrowBandAfterStep(15)
EikoSolv.computeField()
#EikoSolv.fieldValues.show(grid=True, colorscale_name = 'magenta')

fig = go.Figure()
EikoSolv.fieldValues.add3Dplot(fig,color=[255,0,0,0.5])

Opt2=dict(method="FME")

FMEsolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt2)
#EikoSolv.showNarrowBandAfterStep(15)
FMEsolv.computeField()
FMEsolv.fieldValues.add3Dplot(fig,color=[0,255,0,0.5])

MySol.add3Dplot(fig,color=[0,0,255,0.5])

fig.show()

"""
Diff = Mesh2D.VertexValueMap(MyMesh)
Diff.values = np.abs(np.array(EikoSolv.fieldValues.values)/np.array(MySol.values))
print("Max difference : ", max(Diff.values), "Max value of sol : ", max(MySol.values))
print("Total difference : ", sum(Diff.values)/len(Diff.values))
#Diff.show(grid=True, colorscale_name = 'magenta')
PercentDiff = Mesh2D.VertexValueMap(MyMesh)
PercentDiff.values = [Diff.values[i] if MySol.values[i] == 0 else (100*Diff.values[i]/MySol.values[i] if (Diff.values[i]/MySol.values[i] < 0.2) else 0) for i in range(len(MyMesh.vertices)) ]

maxPercent = 0
argMaxPercent = 0
for vertex in range(len(MyMesh.vertices)):
    if PercentDiff[vertex] > maxPercent:
        maxPercent = PercentDiff[vertex]
        argMaxPercent = vertex
print("Max percent difference : ", maxPercent, "At vertex : ", argMaxPercent)
PercentDiff.show(grid=True, colorscale_name = 'magenta')
"""


import hughes2d.Plotter as plth

fig2 = go.Figure()

plth.addVectorFieldPlot(fig2, MyMesh.vertices, MyMesh.triangles, EikoSolv.fieldValues.computeVertexGradientFlow(), color=[255,0,0,0.5])
plth.addVectorFieldPlot(fig2, MyMesh.vertices, MyMesh.triangles, MySol.computeVertexGradientFlow(), color=[0,0,255,0.5])
plth.addVectorFieldPlot(fig2, MyMesh.vertices, MyMesh.triangles, FMEsolv.fieldValues.computeVertexGradientFlow(), color=[0,255,0,0.5])

fig2.show()

#plth.plotVectorField(MyMesh.vertices, MyMesh.triangles, EikoSolv.fieldValues.computeVertexGradientFlow())

"""
plth.plotVectorField(MyMesh.vertices, MyMesh.triangles, MySol.computeGradientFlow())
plth.plotVectorField(MyMesh.vertices, MyMesh.triangles, MySol.computeVertexGradientFlow())
"""
