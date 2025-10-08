from hughes2d.Mesh2D import *

bigSquarePoints = [[0,0],[5,0],[5,5],[0,5]]

smallSquarePoints = [[2,2],[2,3],[3,3],[3,2]]

Exit1 = [[0,1],[0,2]]
Exit2 = [[2.25,3],[2.75,3]]
Exit3 = [[5,1],[5,2]]

Innerwall = [[1,1],[1,4]]

#Domain construction -----------------------------
Domain1 = NonConvexDomain(bigSquarePoints)
Domain1.addWall(smallSquarePoints, cycle=True)
Domain1.addWall(Innerwall, cycle=False)
Domain1.addExit(Exit1)
Domain1.addExit(Exit2)
Domain1.addExit(Exit3)
#Domain1.show()

#Definition of the mesh and generation from the domain -------------
MyMesh = Mesh()
MyMesh.generateMeshFromDomain(Domain1, 0.1)
MyMesh.show()

#Saving the mesh
MyMesh.saveToJson("examples/01-Mesh_generation/outputs/square")

#Loading the mesh
MyMesh.loadFromJson("examples/01-Mesh_generation/outputs/square_mesh.json")
