from hughes2d import *

#Construction du domaine--------------------------------
MyDomain = Mesh2D.NonConvexDomain([[0,0],[0,1],[1,1],[1,0]])
MyDomain.addExits([[[1,0],[1,1]]])
MyDomain.show()

#construction du Mesh--------------------------------------
MyMesh = Mesh2D.Mesh()
MyMesh.generateMeshFromDomain(MyDomain, 0.01)
MyMesh.show()
MyMesh.saveToJson("gettingStartedSimu")


#Construction de la donnée initiale ---------------------------------------
MyMap = Mesh2D.CellValueMap(MyMesh)
MyMap.generateRandom()
MyMap.show()

#Setting the options for the simulation-----------------------------------------
opt = dict(model = "hughes",
            filename = "gettingStartedSimu",
            save = True,
            verbose = True
            )

#Creating the solver and computing---------------------------------------------------
Solver = Splitting.PedestrianSolver(MyMesh, 0.01,0.01, initialDensity = MyMap, options=opt)
Solver.computeUntilEmpty(100)

#Converting the data to a mp4 video------------------------------------------
Plotter.convertToMP4("gettingStartedSimu", limits=[[0,1],[0,1]])
