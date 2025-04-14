from hughes2d import *



MyMesh = Mesh2D.Mesh()

MyMesh.loadFromJson("configZone2_mesh.json")

#MyMesh.show()

MyMap = Mesh2D.CellValueMap(MyMesh)
#MyMap.generateRandom()
for i,t in enumerate(MyMesh.barycenters):
    if(t[0] < 3):
        MyMap.values[i] = 0.99
MyMap.show()

opt = dict(constantDirectionField = False,
            filename = "../../../Visualiseur/data/examples/additional_computations_test3",
            save = True,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.000001},
            eikoSolver = {  'constrained' : False,
                            'NarrowBandDepth' : 2},
            additional_computations = { 'zones_mean_density' : True })

Solver = Splitting.HughesScheme(MyMesh, 0.003,0.0025, initialDensity = MyMap, options=opt)

#for i in range(5):
#Solver.computeStepsAndShow(2)
#Solver.directions[-1].checkGradientValidity()
#Solver.directions[-1].showVectorField()
Solver.computeSteps(5000)
MyMesh.saveToJson(opt['filename'])
#Solver.saveToJson()
