from hughes2d import *


#filename = "mesh_FF.msh"

#MyMesh = Mesh2D.Mesh()

#MyMesh.importMeshFromMshFreeFem(filename)
#MyMesh.saveToJson("config1")

MyMesh = Mesh2D.Mesh()

MyMesh.loadFromJson("config1_mesh.json")

#MyMesh.show()

MyMap = Mesh2D.CellValueMap(MyMesh)
#MyMap.generateRandom()
MyMap.setConstantCircle(center = [3.5,3.5], radius = 2, value=0.6)
MyMap.show()

opt = dict(constantDirectionField = False,
            filename = "../../../Visualiseur/data/Config1_10",
            save = True,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.000001},
            eikoSolver = {  'constrained' : False,
                            'NarrowBandDepth' : 2},
            additional_computations = { 'total_mass' : True })

Solver = Splitting.HughesScheme(MyMesh, 0.03,0.035, initialDensity = MyMap, options=opt)
#for i in range(5):
#Solver.computeStepsAndShow(2)
#Solver.directions[-1].checkGradientValidity()
#Solver.directions[-1].showVectorField()
Solver.computeSteps(1200)
#MyMesh.saveToJson(opt['filename'])
#Solver.saveToJson()
