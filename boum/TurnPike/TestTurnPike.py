from hughes2d import *

import numpy as np

MyMesh = Mesh2D.Mesh()

MyMesh.loadFromJson("testTP_mesh.json")

Xid = [3,6] #np.arange(3,7,1)

for xid in Xid:

    MyMap = Mesh2D.CellValueMap(MyMesh)
    """if(xid == 0):
        MyMap.generateRandom()
    elif(i==3):
        MyMap.setConstantCircle([xid,2.5],2.4,0.7)
    else:
        for i,t in enumerate(MyMesh.barycenters):
            if(1 <= t[0] < 3):
                MyMap.values[i] = 1
            else:
                MyMap.values[i] = 0"""
    MyMap.generateRandom()
    #MyMap.show()


    print("Hughes model number : ", xid)

    opt = dict(model = "hughes",
                filename = "data/TestTPID0_"+str(np.round(xid,1)),
                save = True,
                verbose = True,
                lwrSolver = {   'convexFlux' : True,
                                'anNum' : "dichotomy",
                                'method' : "midVector",
                                'ApproximationThreshold' : 0.000001},
                eikoSolver = {  'constrained' : True,
                                'NarrowBandDepth' : 2},
                additional_computations = { 'total_mass' : True })

    Solver = Splitting.PedestrianSolver(MyMesh, 0.015,0.03, initialDensity = MyMap, costFunction =(lambda x : 1+x), options=opt)

    Solver.computeUntilEmpty(7000)
