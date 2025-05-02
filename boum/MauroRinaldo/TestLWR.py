from hughes2d import *
import hughes2d.Plotter as plth
import random as alea
import plotly

import Solver as MRsolv

import numpy as np


MyMesh = Mesh2D.Mesh()

MyMesh.loadFromJson("configZone_mesh.json")

#MyMesh.show()

MyMap = Mesh2D.CellValueMap(MyMesh)
#MyMap.generateRandom()
MyMap.setConstantCircle([1.5,1],0.7,0.7)
#MyMap.show()

radius = 0.2

epsilon = 0.5

def conv_func(dens, dx, dy):
    return dens*(1-(dx/radius)**2)**3 * (1-(dy/radius)**2)**3

ConstantVectorField = []

for P in MyMesh.barycenters:
    if(1<=P[1]<=2):
        ConstantVectorField.append([1,0])
    elif(P[1]>2):
        V = [3-P[0],2-P[1]]
        ConstantVectorField.append([V[0]/np.sqrt(V[0]**2 + V[1]**2), V[1]/np.sqrt(V[0]**2 + V[1]**2)])
    elif(P[1]<1):
        V = [3-P[0],1-P[1]]
        ConstantVectorField.append([V[0]/np.sqrt(V[0]**2 + V[1]**2), V[1]/np.sqrt(V[0]**2 + V[1]**2)])


opt = dict(constantDirectionField = True,
            filename = "data/Test2champConstant",
            save = True,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001},
            eikoSolver = {  'constrained' : False,
                            'NarrowBandDepth' : 2})

plth.plotVectorField(MyMesh.vertices, MyMesh.triangles, ConstantVectorField)

MySolver = MRsolv.MauroRinaldoScheme(MyMesh, 0.005,0.05, initialDensity = MyMap, constantVectorField=ConstantVectorField, options=opt)

MySolver.computeSteps(2000)

MySolver.lastDensity.show()





#Solver.computeSteps(100)
#Solver.directions[0].show()
#Solver.directions[0].showVectorField()
#Solver.saveToJson()
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
