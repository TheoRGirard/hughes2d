from hughes2d.Mesh2D import *


filename = "examples/06-File_handling/dxf/config_simple.dxf"

MyDomain = NonConvexDomain()
MyDomain.importFromDXF(filename)
MyDomain.show()

MyMesh = Mesh()
MyMesh.generateMeshFromDomain(MyDomain, 10)
MyMesh.show()
