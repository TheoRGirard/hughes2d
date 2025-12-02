import numpy as np
import plotly.graph_objects as go

from hughes2d import (CellValueMap, EikoSolver, Mesh, NonConvexDomain,
                      VertexValueMap)

dx = 0.05

#Construction of the mesh
Exits = [[[12,3],[12,4]]]
ListOuterWall2 = [[0,0],[10,0],[10,3],[12,3],
                  [12,4],[10,4],[10,5],[0,5]]
MyDomain = NonConvexDomain(ListOuterWall2)
MyDomain.add_exits(Exits)
MyMesh = Mesh()
MyMesh.generate_mesh_from_domain(MyDomain,dx)
MyMap = CellValueMap(MyMesh)

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

MySol = VertexValueMap(MyMesh)
MySol.values = ExplicitSolution

Opt=dict(method = "FMT", constrained = False, NarrowBandDepth = 2)

FMTuSolv = EikoSolver(MyMesh, density_map = MyMap, opt = Opt)
FMTuSolv.compute_field()

Opt3=dict(method = "FMT", constrained = True, NarrowBandDepth = 2)

FMTCSolv = EikoSolver(MyMesh, density_map = MyMap, opt = Opt3)
FMTCSolv.compute_field()

Opt2=dict(method = "FME", constrained = True, NarrowBandDepth = 2)

FMEsolv = EikoSolver(MyMesh, density_map = MyMap, opt = Opt2)
FMEsolv.compute_field()


fig = go.Figure()

FMTuSolv.field_values.add_3d_plot(fig,color=[255,127,0,0.5])
FMTCSolv.field_values.add_3d_plot(fig,color=[0,255,0,0.5])
FMEsolv.field_values.add_3d_plot(fig,color=[0,0,255,0.5])

MySol.add_3d_plot(fig,color=[255,0,0,0.5])

fig.show()
