from hughes2d import *

import numpy as np

MyMesh = Mesh2D.Mesh()

MyMesh.loadFromJson("test5_mesh.json")

#MyMesh.show()

MyMap = Mesh2D.CellValueMap(MyMesh)
#MyMap.generateRandom()

MyMap.setConstantCircle([7,2.5],2.4,0.7)
"""for i,t in enumerate(MyMesh.barycenters):
    if(6 <= t[0] < 8):
        MyMap.values[i] = 1"""
MyMap.show()

#raise ValueError("STOP")

print("Hughes model : ")

opt = dict(model = "hughes",
            filename = "data/Test5demo",
            save = True,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.000001},
            eikoSolver = {  'constrained' : True,
                            'NarrowBandDepth' : 2},
            additional_computations = { 'total_mass' : True, 'zones_mean_density' : True, 'max_density' : True })

Solver = Splitting.PedestrianSolver(MyMesh, 0.01,0.02, initialDensity = MyMap, costFunction =(lambda x : 1+5*x), options=opt)

Solver.computeUntilEmpty(7000)



"""

ConstantVectorField = []

for P in MyMesh.barycenters:
    if(P[0]<=5):
        if(3<= P[1]<= 4):
            ConstantVectorField.append([-1,0])
        elif(P[1]>4):
            V = [0-P[0],4-P[1]]
            ConstantVectorField.append([V[0]/np.sqrt(V[0]**2 + V[1]**2), V[1]/np.sqrt(V[0]**2 + V[1]**2)])
        elif(P[1]<3):
            V = [0-P[0],3-P[1]]
            ConstantVectorField.append([V[0]/np.sqrt(V[0]**2 + V[1]**2), V[1]/np.sqrt(V[0]**2 + V[1]**2)])
    else:
        if(3<= P[1]<= 4):
            ConstantVectorField.append([1,0])
        elif(P[1]>4):
            V = [10-P[0],4-P[1]]
            ConstantVectorField.append([V[0]/np.sqrt(V[0]**2 + V[1]**2), V[1]/np.sqrt(V[0]**2 + V[1]**2)])
        elif(P[1]<3):
            V = [10-P[0],3-P[1]]
            ConstantVectorField.append([V[0]/np.sqrt(V[0]**2 + V[1]**2), V[1]/np.sqrt(V[0]**2 + V[1]**2)])


print("LWR pur model : ")

opt = dict(model = "constantDirectionField",
            filename = "data/LWRconstantFieldID2",
            save = True,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.000001},
            eikoSolver = {  'constrained' : False,
                            'NarrowBandDepth' : 2},
            additional_computations = { 'total_mass' : True, 'zones_mean_density' : True, 'max_density' : True })

Solver = Splitting.PedestrianSolver(MyMesh, 0.01,0.02, initialDensity = MyMap, directions = ConstantVectorField, options=opt)

Solver.computeSteps(4000)



print("Colombo Garavello : ")

opt = dict(model = "colombo-garavello",
            CGparameters = {"radius" : 0.3, "epsilon" : 0.5},
            filename = "data/ColomboGaravelloID2",
            save = True,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.000001},
            eikoSolver = {  'constrained' : False,
                            'NarrowBandDepth' : 2},
            additional_computations = { 'total_mass' : True, 'zones_mean_density' : True, 'max_density' : True })

Solver = Splitting.PedestrianSolver(MyMesh, 0.01,0.02, initialDensity = MyMap, directions = ConstantVectorField, options=opt)

Solver.computeSteps(4000)
"""
