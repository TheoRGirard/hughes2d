from hughes2d import *
import random as alea
import plotly

import numpy as np


Exits = [[[5,0],[5,2]]]

ListOuterWall2 = [[0,0],[5,0],[5,0.7],[5,1.3],[5,2],[0,2]]


MyDomain = Mesh2D.NonConvexDomain(ListOuterWall2)
MyDomain.addExits(Exits)

"""r = 0.2
N = 10
InnerWall = [[4.6+ r*np.cos(2*i/N*np.pi), 1+ r*np.sin(2*i/N*np.pi)] for i in range(N)]

MyDomain.addWall(InnerWall,cycle=True)"""

MyDomain.show()

MyMesh = Mesh2D.Mesh(MyDomain, 0.005)

MyMesh.triangulate()
#MyMesh.show()

MyMap = Mesh2D.CellValueMap(MyMesh)
#MyMap.generateRandom()
for i in range(len(MyMap)):
    if(MyMesh.meshDict["vertices"][MyMesh.meshDict["triangles"][i][0]][0] < 1.5 and
        MyMesh.meshDict["vertices"][MyMesh.meshDict["triangles"][i][0]][0] > 0.2 and
        MyMesh.meshDict["vertices"][MyMesh.meshDict["triangles"][i][0]][1] < 1.6 and
        MyMesh.meshDict["vertices"][MyMesh.meshDict["triangles"][i][0]][1] > 0.4):
        MyMap[i] = 0.5
    else:
        MyMap[i] = 0
MyMap.show()
#print(MyMesh.meshDict)
#print(len(MyMesh.meshDict['segments']))
#print(len(MyMesh.meshDict['vertices']))

"""

MyEmptyMap = CellValueMap(MyMesh)

Opt=dict(constrained = False, NarrowBandDepth = 2)

MyEikoSolv = EikoSolver(MyMesh, DensityMap=MyEmptyMap, opt=Opt)
MyEikoSolv.computeField()
MyEikoSolv.fieldValues.showVectorField()

MyVectors = MyEikoSolv.fieldValues.computeGradientFlow()"""

MyVectors = [[1,0] for i in MyMesh.meshDict['triangles']]

opt = dict(constantDirectionField = False,
            filename = "data/HughesCouloirTmap",
            save = True,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "tmap",
                            'ApproximationThreshold' : 0.0001},
            eikoSolver = {  'constrained' : False,
                            'NarrowBandDepth' : 2})

Solver = Splitting.HughesScheme(MyMesh, 0.0025,0.005, initialDensity = MyMap, directions = MyVectors, options=opt)
#for i in range(20):
    #Solver.computeStepsAndShow(30)

#Solver.computeStepsAndShow(20)

Solver.computeSteps(100)
#Solver.directions[0].show()
#Solver.directions[0].showVectorField()
Solver.saveToJson()
#Solver.computeStepsAndShow(1)
#Solver.show()
#Solver.computeStepsAndShow(5)

"""print(MyMesh.meshDict['triangles'][MyMesh.ListOfPairTriangles[0][0]])
print(MyMesh.meshDict['triangles'][MyMesh.ListOfPairTriangles[0][1]])
print(MyMesh.TriangleWithEdgeList[MyMesh.ListOfPairTriangles[0][0]])
print(MyMesh.TriangleWithEdgeList[MyMesh.ListOfPairTriangles[0][1]])
print(MyMesh.outerNormalVectByTriangles[MyMesh.ListOfPairTriangles[2][0]])
print(MyMesh.outerNormalVectByTriangles[MyMesh.ListOfPairTriangles[2][1]])"""

#print(MyMesh.getExitVertices())

#print(ComputeOuterNormalUnitVect([5.97,-10.47],[6,-10.75],[5.72,-10.33]))
#MyDomain2 = NonConvexDomain([[5.97,-10.47],[6,-10.75],[5.72,-10.33]])
#MyDomain2.show()

#Solver = EikoSolver(MyMesh, exitVertices=MyMesh.getExitVertices())
#Solver.computeField()
#print(Solver.fieldValues)
#Solver.fieldValues.show()
#Solver.fieldValues.showVectorField()
