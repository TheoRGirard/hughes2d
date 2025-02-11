from hughes2d import *

filename = "mesh_FF.msh"

MyMesh = Mesh2D.Mesh()

MyMesh.loadFromJson("config1_mesh.json")

MyMesh.show()
