from hughes2d.Mesh2D import *

filename = "examples/06-File_handling/FreeFem/mesh_FF.msh"

MyMesh = Mesh()
MyMesh.importMeshFromMshFreeFem(filename)
MyMesh.show()
