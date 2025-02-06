from hughes2d.Mesh2D import *

filename = "TestLibreCad.dxf"

MyDomain = NonConvexDomain()

MyDomain.importFromDXF(filename)

MyDomain.show()


MyMesh = Mesh()

MyMesh.generateMeshFromDomain(MyDomain, 15)

MyMesh.show()
