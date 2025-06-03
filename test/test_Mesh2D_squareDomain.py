from hughes2d.Mesh2D import *
import pytest


#Tests for the NonConvexDomain object -------------------------------------------

bigSquarePoints = [[0,0],[5,0],[5,5],[0,5]]

smallSquarePoints = [[2,2],[2,3],[3,3],[3,2]]

Exit1 = [[0,1],[0,2]]
Exit2 = [[2.25,3],[2.75,3]]

InWall = [[1,1],[3,3]]

Domain1 = NonConvexDomain(bigSquarePoints)

@pytest.mark.parametrize("point,edge,triangle", [([0,0],[[-1,-1],[1,1]],[[-1,-1],[1,0],[0,1]])])
def test_NCD_staticBelongs(point, edge, triangle):
    assert belongTriangle(point,triangle)
    assert NonConvexDomain.belongSegment(point,edge)

def test_NCD_contains():
    assert [2,2] in Domain1
    assert [-1,2] not in Domain1
    assert [6,6] not in Domain1

def test_NCD_hasBoundaryPoint():
    assert Domain1.hasBoundaryPoint([0,3])
    assert not Domain1.hasBoundaryPoint([1,1])

def test_NCD_addWall():
    with pytest.raises(ValueError):
        Domain1.addWall([[0,0],[1,-1]], cycle=False)

Domain1.addWall(smallSquarePoints, cycle=True)

def test_NCD_holes():
    assert [2.5,2.5] in Domain1
    assert not [6,0] in Domain1
    #Caution : a hole in the domain does NOT exclude points in the __contains__ method at the moment.

def test_NCD_addBoundaryPoint():
    with pytest.raises(ValueError):
        Domain1.addBoundaryPoint([2.5,2.5])
        Domain1.addBoundaryPoint([2.5,3])
Domain1.addBoundaryPoint([2,5])

def test_NCD_addWallPoint():
    with pytest.raises(ValueError):
        Domain1.addWallPoint([2.5,2.5])
        Domain1.addWallPoint([0,2.5])
Domain1.addWallPoint([2.5,3])

def test_NCD_hasWallPoint():
    assert Domain1.hasWallPoint([2.5,3])
    assert not Domain1.hasWallPoint([0,3])
    #Caution : wall stands for inner wall. For the outer boundary, everything that is not considered as an exit will be treated as a wall.
    assert not Domain1.hasWallPoint([1,1])

Domain1.addExits([Exit1, Exit2])

def test_NCD_hasExitPoint():
    assert Domain1.hasExitPoint([2.5,3])
    assert not Domain1.hasExitPoint([2.5,2])
    assert not Domain1.hasExitPoint([1,1])

def test_NCD_zones():
    Domain1.addZone("Bottom side", [[0,0],[0,2],[5,2],[5,0]])
    with pytest.raises(NameError):
        Domain1.addZone("Bottom side", [[0,0],[0,2],[5,2],[5,0]])


@pytest.mark.parametrize("edge,result", [([[2.5,3],[2.6,3]],True),([[2.5,2],[2.6,2]],True),([[2.5,0],[2.6,0]],False)])
def test_NCD_hasWallEdge(edge, result):
    assert Domain1.hasWallEdge(edge) == result

@pytest.mark.parametrize("edge,result", [([[2.5,3],[2.6,3]],True),([[2.5,2],[2.6,2]],False),([[2.5,0],[2.6,0]],False)])
def test_NCD_hasExitEdge(edge, result):
    assert Domain1.hasExitEdge(edge) == result

@pytest.mark.parametrize("edge,result", [([[2.5,3],[2.6,3]],False),([[2.5,2],[2.6,2]],False),([[2.5,0],[2.6,0]],True)])
def test_NCD_hasOuterEdge(edge, result):
    assert Domain1.hasOuterEdge(edge) == result



#Tests for the Mesh object ----------------------------------------------------

Mesh1 = Mesh()
Mesh1.generateMeshFromDomain(Domain1,0.1)

def test_Mesh_zones():
    with pytest.raises(NameError):
        Mesh1.inZone([1,1], "Left side")

    Mesh1.addConvexZone("Left side", [[0,0],[2,0],[2,5],[0,5]])
    with pytest.raises(NameError):
        Mesh1.addConvexZone("Left side", [[0,0],[2,0],[2,5],[0,5]])

    assert Mesh1.inZone([1,1], "Left side")

    Mesh1.computeZonesTriangles()
    assert len(Mesh1.zones["Left side"]['triangles']) > 0

def test_generateMesh():
    assert max(Mesh1.cellAreas) <  0.1

def test_saveLoadJson():
    Mesh1.saveToJson("test/testMesh1")
    Mesh2 = Mesh()
    Mesh2.loadFromJson("test/testMesh1_mesh.json")

    for i in range(len(Mesh1.triangles)):
        for j in range(3):
            assert Mesh1.triangles[i][j] == Mesh2.triangles[i][j]

def test_exits():
    for indexP in Mesh1.exitVertices:
        assert Domain1.hasExitPoint(Mesh1.vertices[indexP])
    for edge in Mesh1.exitEdges:
        exit = [    [Mesh1.vertices[Mesh1.edges[edge][0]][0], Mesh1.vertices[Mesh1.edges[edge][0]][1]],
                    [Mesh1.vertices[Mesh1.edges[edge][1]][0], Mesh1.vertices[Mesh1.edges[edge][1]][1]] ]
        assert Domain1.hasExitEdge(exit)

def test_walls():
    for indexP in Mesh1.boundaryPoints:
        assert Domain1.hasBoundaryPoint(Mesh1.vertices[indexP]) or Domain1.hasWallPoint(Mesh1.vertices[indexP])
    for edge in Mesh1.wallEdges:
        wall = [    [Mesh1.vertices[Mesh1.edges[edge][0]][0], Mesh1.vertices[Mesh1.edges[edge][0]][1]],
                    [Mesh1.vertices[Mesh1.edges[edge][1]][0], Mesh1.vertices[Mesh1.edges[edge][1]][1]] ]
        assert Domain1.hasWallEdge(wall) or Domain1.hasOuterEdge(wall)



#Tests for the CellValueMap object--------------------------------------------
Map1 = CellValueMap(Mesh1)

for i in range(len(Map1)):
    Map1[i] = 0.1

def test_algebra_cellvaluemap():
    Map2 = CellValueMap(Mesh1)
    for i in range(len(Map2)):
        Map2[i] = 0.3
    Map3 = Map1*2 + Map2
    for i in range(len(Map3)):
        assert Map3[i] == 0.5
    Map3 += 0.25
    print(type(Map3))
    for i in range(len(Map3)):
        assert Map3[i] == 0.75

def test_sets_cellvaluemap():
    Map2 = CellValueMap(Mesh1)
    for i in range(len(Map2)):
        assert Map2[i] == 0
    Map2.generateRandom(variability = 0.4,mean = 0.3)
    for i in range(len(Map2)):
        assert Map2[i] > 0

    Map1.setConstantCircle(center=[1,1], radius = 1, value = 0.5)
    for i in range(len(Map1)):
        if np.sqrt( (Mesh1.barycenters[i][0]-1)**2 + (Mesh1.barycenters[i][1]-1)**2 ) < 1:
            assert Map1[i] == 0.5
        else:
            assert Map1[i] == 0.1

def test_convolutionOverSquareBall():
    Mesh2 = Mesh()
    Mesh2.generateMeshFromDomain(Domain1, 0.02,da=10)
    Map = CellValueMap(Mesh2)
    for i in range(len(Map)):
        Map[i] = 2
    Result = Map.convolutionOverSquareBall(radius=0.2, conv_func = (lambda x,y,z:2))
    for i in range(len(Mesh2.vertices)):
        if (Mesh2.vertices[i][0] > 0.2 and Mesh2.vertices[i][0] < 4.8
                and Mesh2.vertices[i][1] > 0.2 and Mesh2.vertices[i][1] < 4.8
                and (Mesh2.vertices[i][0] < 1.8 or Mesh2.vertices[i][0] > 3.2)
                and (Mesh2.vertices[i][1] < 1.8 or Mesh2.vertices[i][1] > 3.2)) :
            assert np.abs(Result[i] - 2*0.4**2) < float(1e-1)
        else:
            assert Result[i] < 2*0.4**2 + float(1e-1)
