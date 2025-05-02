from hughes2d import *

#filename = "mesh_FF.msh"
#filename = "TestExportSansTrou_FF.msh"

filename = "config 1.dxf"

MyDomain = Mesh2D.NonConvexDomain()

MyDomain.importFromDXF(filename)

MyDomain.show()


MyMesh = Mesh2D.Mesh()

MyMesh.generateMeshFromDomain(MyDomain, 0.5)

MyMesh.show()

MyMesh.saveToJson("configTest")

MyMesh.exportMeshMshFreeFem("configTest.msh")
