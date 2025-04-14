from hughes2d import *

import numpy as np

#filename = "mesh_FF.msh"
#filename = "TestExportSansTrou_FF.msh"


MyDomain = Mesh2D.NonConvexDomain()

MyDomain = Mesh2D.NonConvexDomain([[0,0],[3,0],[3,1],[5,1],[5,2],[3,2],[3,3],[0,3]])

MyDomain.addExit([[5,1],[5,2]])

#MyDomain.show()

for dx in np.arange(0,5,0.025):
    MyDomain.addZone("column_start_"+str(round(dx,3)), [[round(dx,3),0],[round(dx+0.1,3),0],[round(dx+0.1,3),3],[round(dx,3),3]])


MyMesh = Mesh2D.Mesh()

MyMesh.generateMeshFromDomain(MyDomain, 0.0025)

#print(MyMesh.zones)

MyMesh.show()

MyMesh.saveToJson("configZone2")
