from hughes2d.Mesh2D import *
import pytest

bigSquarePoints = [[0,0],[5,0],[5,5],[0,5]]

smallSquarePoints = [[2,2],[2,3],[3,3],[3,2]]

Exit1 = [[0,1],[0,2]]
Exit2 = [[2.25,3],[2.75,3]]

InWall = [[1,1],[3,3]]

Domain1 = NonConvexDomain(bigSquarePoints)

#Domain1.addWall(smallSquarePoints, cycle=True)

Domain1.addExit(Exit1)
#Domain1.addExit(Exit2)

#Domain1.show()

Mesh1 = Mesh()

Mesh1.generateMeshFromDomain(Domain1, 0.1)

Mesh1.show()

#print(Mesh1.exitEdges, Mesh1.wallEdges)

Mesh1.exportMeshMshFreeFem("TestExportSansTrou_FF.msh")
