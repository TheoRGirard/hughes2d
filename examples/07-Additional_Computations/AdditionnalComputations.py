from pathlib import Path

from hughes2d import CellValueMap, Mesh, NonConvexDomain, PedestrianSolver

MyMesh = Mesh()

MyDomain = NonConvexDomain([[0,0],[2,0],[2,1],[0,1]])
MyDomain.add_exit([[2,0],[2,1]])
MyDomain.add_zone("Corridor_left", [[0,0],[0,1],[1,1],[1,0]])
MyDomain.add_zone("Corridor_right", [[2,0],[2,1],[1,1],[1,0]])
MyMesh.generate_mesh_from_domain(MyDomain, 0.01)

MyMap = CellValueMap(MyMesh)
for i,B in enumerate(MyMesh.barycenters):
    if(B[0] < 0.4):
        MyMap.values[i] = 0.5


opt = dict(constantDirectionField = False,
            filename = Path(__file__).parent / "data" / "output",
            save = True,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001},
            eikoSolver = {  'constrained' : False,
                            'NarrowBandDepth' : 2},
            additional_computations = { 'zones_mean_density' : True, 'total_mass' : True })

Solver = PedestrianSolver(MyMesh, 0.003, initial_density = MyMap, options=opt)

Solver.compute_until_empty(5000)
Solver.save_to_json()
