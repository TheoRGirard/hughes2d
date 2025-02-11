from hughes2d import *

filename = "mesh_FF.msh"
#filename = "TestExportSansTrou_FF.msh"

MyMesh = Mesh2D.Mesh()

MyMesh.importMeshFromMshFreeFem(filename)

MyMesh.saveToJson("config1")
