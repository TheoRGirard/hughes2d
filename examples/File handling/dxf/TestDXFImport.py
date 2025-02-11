from hughes2d.Mesh2D import *

#filename = "TestLibreCad.dxf"
filename = "config 1 pillar.dxf"

MyDomain = NonConvexDomain()

MyDomain.importFromDXF(filename)

MyDomain.show()


MyMesh = Mesh()

MyMesh.generateMeshFromDomain(MyDomain, 10)

MyMesh.show()
