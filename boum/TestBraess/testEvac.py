from hughes2d import *

import numpy as np

#filename = "mesh_FF.msh"
#filename = "TestExportSansTrou_FF.msh"

disc = 19
r = 0.7
"""for x in np.arange(3,3.6,0.1):

    print("Poteau à ", x)"""

x = 3.5

MyDomain = Mesh2D.NonConvexDomain()

MyDomain = Mesh2D.NonConvexDomain([[0,0],[5,0],[5,2],[7,2],[7,3],[5,3],[5,5],[0,5]])

MyDomain.addExit([[7,2],[7,3]])

#MyDomain.show()


MyDomain.addZone("room", [[0,0],[5,0],[5,5],[0,7]])
MyDomain.addZone("corridor", [[5,2],[7,2],[7,3],[5,3]])



MyDomain.addWall([[x+r*np.cos(2*np.pi*i/(disc)),2.5+r*np.sin(2*np.pi*i/(disc))] for i in range(disc)], cycle=True)

MyDomain.show()


MyMesh = Mesh2D.Mesh()

MyMesh.generateMeshFromDomain(MyDomain, 0.01)

#print(MyMesh.zones)
MyMesh.saveToJson("data/HughesX3.5")
MyMesh.show()


MyMap = Mesh2D.CellValueMap(MyMesh)
#MyMap.generateRandom()

#print([ index for index,Pair in enumerate(MyMesh.pairsOfTriangles) if (len(Pair) == 1 and index not in MyMesh.exitEdges and index not in My Mesh.wallEdges)])
#print(len(MyMesh.exitEdges) + len(MyMesh.wallEdges))

#MyMap.setConstantCircle([3.5,3.5],2,0.7)
for i,t in enumerate(MyMesh.barycenters):
    if(t[0] < 2 and (1.5<t[1]< 3.5)):
        MyMap.values[i] = 1
MyMap.show()

print("Hughes model : ")

opt = dict(model = "hughes",
            filename = "data/HughesX3.5",
            save = True,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.000001},
            eikoSolver = {  'method' : 'FMT',
                            'constrained' : True,
                            'NarrowBandDepth' : 2}
                            )

Solver = Splitting.PedestrianSolver(MyMesh, 0.01,0.01, initialDensity = MyMap, costFunction =(lambda x : 1+x), options=opt)

Solver.computeUntilEmpty(7000)
