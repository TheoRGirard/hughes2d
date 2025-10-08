from hughes2d.Mesh2D import *

#From : https://github.com/nschloe/meshio/blob/main/tests/meshes/msh/insulated-2.2.msh
filename = "examples/06-File_handling/FreeFem/insulated-2.2.msh"


MyMesh = Mesh()
MyMesh.importMeshFromMsh(filename)
MyMesh.show()
