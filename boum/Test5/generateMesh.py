from hughes2d import *

import numpy as np

#filename = "mesh_FF.msh"
#filename = "TestExportSansTrou_FF.msh"

MyDomain = Mesh2D.NonConvexDomain()

MyDomain = Mesh2D.NonConvexDomain([[0,0],[10,0],[10,3],[12,3],[12,4],[10,4],[10,7],[0,7],[0,4],[-2,4],[-2,3],[0,3]])

MyDomain.addExit([[12,3],[12,4]])
MyDomain.addExit([[-2,3],[-2,4]])

#MyDomain.show()


MyDomain.addZone("roomL", [[0,0],[5,0],[5,7],[0,7]])
MyDomain.addZone("roomR", [[10,0],[5,0],[5,7],[10,7]])
MyDomain.addZone("corridorL", [[-2,3],[0,3],[0,4],[-2,4]])
MyDomain.addZone("corridorR", [[10,3],[12,3],[12,4],[10,4]])


#MyDomain.addWall([[x+(0.9-r)+r*np.cos(2*np.pi*i/(disc)),3.5+r*np.sin(2*np.pi*i/(disc))] for i in range(disc)], cycle=True)

MyDomain.show()


MyMesh = Mesh2D.Mesh()

MyMesh.generateMeshFromDomain(MyDomain, 0.03)

#print(MyMesh.zones)

MyMesh.show()


#MyMesh.saveToJson("meshes/BraessX"+str(np.round(x+(0.9-r),2))+"R"+str(np.round(r,2)))
MyMesh.saveToJson("test5")
