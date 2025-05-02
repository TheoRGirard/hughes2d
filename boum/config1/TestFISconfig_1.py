from hughes2d import *

filename = "mesh_FF.msh"

MyMesh = Mesh2D.Mesh()

MyMesh.loadFromJson("config1_mesh.json")

MyMap = Mesh2D.CellValueMap(MyMesh)

MyMap.setConstantCircle(center = [3.5,3.5], radius = 3, value=0.65)
MyMap.show()

speeds = [0.5 + 0.1*i for i in range(15)]

for v in speeds:
    opt = dict(constantDirectionField = False,
                filename = "../../../Visualiseur/data/Config1_FIS_"+str(v),
                save = True,
                verbose = True,
                lwrSolver = {   'convexFlux' : True,
                                'anNum' : "dichotomy",
                                'method' : "midVector",
                                'ApproximationThreshold' : 0.0001},
                eikoSolver = {  'constrained' : False,
                                'NarrowBandDepth' : 2})

    print("In file : ", opt['filename'])
    speed = lambda x : v*(1-x)

    Solver = Splitting.HughesScheme(MyMesh, 0.025,0.035, initialDensity = MyMap, speedFunction = speed, options=opt)
    Solver.computeSteps(1000)
