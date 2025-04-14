from hughes2d import *

import numpy as np

MyMesh = Mesh2D.Mesh()

MyMesh.loadFromJson("configZone_mesh.json")

#MyMesh.show()

MyMap = Mesh2D.CellValueMap(MyMesh)
#MyMap.generateRandom()
MyMap.setConstantCircle([1.5,1],0.7,0.7)
#MyMap.show()

radius = 0.2

def conv_func(dens, dx, dy):
    return dens*(1-(dx/radius)**2)**3 * (1-(dy/radius)**2)**3

MyConv = Mesh2D.VertexValueMap(MyMesh)
MyConv.values = MyMap.convolutionOverSquareBall(radius, conv_func)

#MyConv.show()
epsilon = 0.5
MyConv.showVectorField(normalize=True, normalization = (lambda x,y: radius**2 * (np.sqrt(1 + x**2 + y**2))/epsilon) )
