from hughes2d.Mesh2D import *

#filename = "TestLibreCad.dxf"
filename = "config 2 constraint.dxf"

MyDomain = NonConvexDomain()

MyDomain.importFromDXF(filename)

#MyDomain.show()

print("Mesh creation")

MyMesh = Mesh()

print("Mesh generation")

MyMesh.generateMeshFromDomain(MyDomain, 0.5)

print("Mesh plot")

MyMesh.show()
