from hughes2d.Mesh2D import *
import pytest

bigSquarePoints = [[0,0],[5,0],[5,5],[0,5]]

smallSquarePoints = [[2,2],[2,3],[3,3],[3,2]]

Exit1 = [[0,1],[0,2]]
Exit2 = [[2.25,3],[2.75,3]]

InWall = [[1,1],[3,3]]

Domain1 = NonConvexDomain(bigSquarePoints)

@pytest.mark.parametrize("point,edge,triangle", [([0,0],[[-1,-1],[1,1]],[[-1,-1],[1,0],[0,1]])])
def test_NCD_staticBelongs(point, edge, triangle):
    assert NonConvexDomain.belongTriangle(point,triangle)
    assert NonConvexDomain.belongSegment(point,edge)

def test_NCD_contains():
    assert [2,2] in Domain1
    assert [-1,2] not in Domain1
    assert [6,6] not in Domain1

def test_NCD_hasBoundaryPoint():
    assert Domain1.hasBoundaryPoint([0,3]) >= 0
    assert Domain1.hasBoundaryPoint([1,1]) == -1

Domain1.addWall(smallSquarePoints, cycle=True)

def test_NCD_holes():
    assert [2.5,2.5] in Domain1
    #Caution : a hole in the domain does NOT exclude points in the __contains__ method at the moment.

def test_NCD_hasWallPoint():
    assert Domain1.hasWallPoint([2.5,3]) >= 0
    assert Domain1.hasWallPoint([0,3]) == -1
    #Caution : wall stands for inner wall. For the outer boundary, everything that is not considered as an exit will be treated as a wall.
    assert Domain1.hasWallPoint([1,1]) == -1

Domain1.addExit(Exit1)
Domain1.addExit(Exit2)
Domain1.show()



#assert Domain1.addWall(InWall) ValueError
