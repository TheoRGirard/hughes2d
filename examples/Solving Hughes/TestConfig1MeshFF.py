from hughes2d import *

#From : https://github.com/nschloe/meshio/blob/main/tests/meshes/msh/insulated-2.2.msh
#filename = "insulated-2.2.msh"
#filename = "TestExportSansTrou_FF.msh"
filename = "mesh_FF.msh"

MyMesh = Mesh2D.Mesh()

MyMesh.importMeshFromMshFreeFem(filename)

#MyMesh.show()

MyMap = Mesh2D.CellValueMap(MyMesh)
#MyMap.generateRandom()
MyMap.setConstantCircle(center = [1,1], radius = 0.5, value=0.5)
MyMap.show()

opt = dict(constantDirectionField = False,
            filename = "../Visualiseur/data/Config1_2",
            save = True,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001},
            eikoSolver = {  'constrained' : False,
                            'NarrowBandDepth' : 2})

Solver = Splitting.HughesScheme(MyMesh, 0.025,0.035, initialDensity = MyMap, options=opt)
#for i in range(5):
#Solver.computeStepsAndShow(2)
#Solver.directions[-1].checkGradientValidity()
#Solver.directions[-1].showVectorField()
Solver.computeSteps(2000)
Solver.saveToJson()
