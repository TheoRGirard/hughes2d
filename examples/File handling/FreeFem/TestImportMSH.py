from hughes2d.Mesh2D import *

#From : https://github.com/nschloe/meshio/blob/main/tests/meshes/msh/insulated-2.2.msh
#filename = "insulated-2.2.msh"
#filename = "TestExport.msh"
filename = "mesh_FF.msh"

MyMesh = Mesh()

MyMesh.importMeshFromMshFreeFem(filename)

MyMesh.show()
