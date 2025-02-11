from Hughes2D import *
import random as alea
import plotly



Exits = [[[8,2],[8,3]]]

ListOuterWall2 = [[0,0],[0,2],[0,3],[0,5],[8,5],[8,0]]



#Construction du domaine--------------------------------
MyDomain = Mesh2D.NonConvexDomain(ListOuterWall2)
MyDomain.addExits(Exits)
#MyDomain.show()

#construction du Mesh--------------------------------------
MyMesh = Mesh2D.Mesh(MyDomain, 0.01)
MyMesh.triangulate()
MyMesh.show()


MyMap = Mesh2D.CellValueMap(MyMesh)
"""for n,T in enumerate(MyMesh.meshDict['triangles']):
    A = MyMesh.meshDict['vertices'][T[0]]
    B = MyMesh.meshDict['vertices'][T[1]]
    C = MyMesh.meshDict['vertices'][T[2]]
    if((A[1]-2.5)*(B[1]-2.5) < 0 or (B[1]-2.5)*(C[1]-2.5) < 0 or (A[1]-2.5)*(C[1]-2.5) < 0):
        MyMap[n] = 0
    else:
        MyMap[n] = 1"""
#MyMap.show()

Opt=dict(constrained = True, NarrowBandDepth = 1)

EikoSolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt)
#EikoSolv.showNarrowBandAfterStep(15)
EikoSolv.computeField()
EikoSolv.fieldValues.show(grid=True, colorscale_name = 'magenta')
"""EikoSolv.fieldValues.showVectorField()
EikoSolv.fieldValues.show3D()"""
