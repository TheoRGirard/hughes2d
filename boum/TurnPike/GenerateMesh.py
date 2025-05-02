from hughes2d import *

import numpy as np

MyDomain = Mesh2D.NonConvexDomain()

MyDomain = Mesh2D.NonConvexDomain([[0,0],[10,0],[10,3],[12,3],[12,4],[10,4],[10,7],[0,7],[0,4],[0,3]])

MyDomain.addExit([[12,3],[12,4]])


MyDomain.show()

MyMesh = Mesh2D.Mesh()

MyMesh.generateMeshFromDomain(MyDomain, 0.03)


MyMesh.show()

MyMesh.saveToJson("testTP")
