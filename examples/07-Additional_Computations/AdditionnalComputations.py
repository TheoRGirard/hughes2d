from hughes2d import *

MyMesh = Mesh2D.Mesh()

MyDomain = Mesh2D.NonConvexDomain([[0,0],[2,0],[2,1],[0,1]])
MyDomain.addExit([[2,0],[2,1]])
MyDomain.addZone("Corridor_left", [[0,0],[0,1],[1,1],[1,0]])
MyDomain.addZone("Corridor_right", [[2,0],[2,1],[1,1],[1,0]])
MyMesh.generateMeshFromDomain(MyDomain, 0.01)

#MyMesh.show()

MyMap = Mesh2D.CellValueMap(MyMesh)
for i,B in enumerate(MyMesh.barycenters):
    if(B[0] < 0.4):
        MyMap.values[i] = 0.5


opt = dict(constantDirectionField = False,
            filename = "examples/07-Additional_Computations/data/output",
            save = True,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001},
            eikoSolver = {  'constrained' : False,
                            'NarrowBandDepth' : 2},
            additional_computations = { 'zones_mean_density' : True, 'total_mass' : True })

Solver = Splitting.PedestrianSolver(MyMesh, 0.003, initialDensity = MyMap, options=opt)

Solver.computeUntilEmpty(5000)
Solver.saveToJson()
