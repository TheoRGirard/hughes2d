"""
2D Mesh utils for triangular meshes on non-convex domains

Last update : 29/08/25

Girard Théo
mail : theo.girard@univ-tours.fr
"""

import numpy as np

from copy import copy, deepcopy
import random as alea

from numpy.typing import NDArray, ArrayLike
from typing import List,Tuple

PointType = List[float]


import triangle as tr

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

try:
    import plotly.figure_factory as ff
    import plotly.graph_objects as go
except ImportError:
    go, ff = None, None

try:
    import ezdxf
except ImportError:
    ezdxf = None

try:
    import meshio
except ImportError:
    meshio = None

try:
    import json
except ImportError:
    json = None


class NonConvexDomain(object):
    """Domain object, it contains :
    - The list of the outer vertices
    - The list of the outer eges
    - The lists of edges and vertices corresponding to the holes inside the domain
    - The list of exit edges
    methods :
    - show : plot the domain
    - addWall/addExit : domain editing tools
    - contains/hasBoundaryPoint : boolean tests for points inside the domain (resp. boundary)
    """

    def __init__(self, outerVerticesList: List[PointType]=[[0,0],[1,0],[0,1]] ) -> None:
        #print([vertexList[edge[0]] for edge in convexHull ])
        self.listOuterVertex = outerVerticesList
        self.outerBoundary = []
        for i in range(len(outerVerticesList)):
            self.outerBoundary.append([i,(i+1)%(len(outerVerticesList))])
        #print(self.outerBoundary)
        self.wallVertices = []
        self.wallEdges = []
        self.wallHolesPoint = []
        self.exitList = []

    def importFromDXF(self, filename: str) -> int:
        """
        Import from a .dxf file. The dxf file can contain three layers :
        - a "domain" layer (mandatory) that contain the outer boundary
        - an "innerWalls" layer that contains all the inner structures of the domain.
        - an "exits" layer that contains all the exits of the domain.

        The package ezdxf must be installed in order to use this method.
        """
        if(ezdxf):
            try:
                doc = ezdxf.readfile(filename)
            except IOError:
                print(f"Not a DXF file or a generic I/O error.")
                return 0
            except ezdxf.DXFStructureError:
                print(f"Invalid or corrupted DXF file.")
                return 0

            print("Extracting data from %s" % filename)
            # helper function
            #def print_entity(e):
            #    print("LINE on layer: %s\n" % e.dxf.layer)
            #    print("start point: %s\n" % e.dxf.start)
            #    print("end point: %s\n" % e.dxf.end)

            self.listOuterVertex = []
            self.outerBoundary = []
            self.wallVertices = []
            self.wallEdges = []
            self.wallHolesPoint = []

            # iterate over all entities in modelspace
            msp = doc.modelspace()

            # entity query for all LINE entities in modelspace
            for e in msp.query("LINE"):
                if(e.dxf.layer == "0" or e.dxf.layer == "domain"):
                    P = [e.dxf.start[0],e.dxf.start[1]]
                    if P not in self.listOuterVertex:
                        numP = len(self.listOuterVertex)
                        self.listOuterVertex.append(P)
                    else:
                        numP = self.addBoundaryPoint(P)

                    Q = [e.dxf.end[0],e.dxf.end[1]]
                    if Q not in self.listOuterVertex:
                        numQ = len(self.listOuterVertex)
                        self.listOuterVertex.append(Q)
                    else:
                        numQ = self.addBoundaryPoint(Q)

                    self.outerBoundary.append([numP, numQ])

            if([len(self.listOuterVertex)-1, 0] not in self.outerBoundary ):
                raise ValueError("Corrupted domain : domain not closed.")
            wallVert = []
            for e in msp.query("LINE"):
                if(e.dxf.layer == "innerWalls"):
                    cycling = 0
                    P = [e.dxf.start[0],e.dxf.start[1]]
                    if P not in wallVert:
                        wallVert.append(P)
                    else:
                        cycling += 1

                    Q = [e.dxf.end[0],e.dxf.end[1]]
                    if Q not in wallVert:
                        wallVert.append(Q)
                    else:
                        cycling += 1

                    if(cycling == 0):
                        #New cycle !
                        wallVert.remove(P)
                        wallVert.remove(Q)
                        if(len(wallVert) > 0):
                            self.addWall(wallVert, cycle =False)
                        wallVert = [P,Q]
                    elif(cycling == 2):
                        self.addWall(wallVert, cycle =True)
                        wallVert = []

            if(len(wallVert) > 0):
                self.addWall(wallVert, cycle =False)

            for e in msp.query("LINE"):
                if(e.dxf.layer == "exits"):
                    P = [e.dxf.start[0],e.dxf.start[1]]
                    Q = [e.dxf.end[0],e.dxf.end[1]]
                    self.addExit([P,Q])

        else:
            raise ImportError("ezdxf not or wrongly installed.")


    def __contains__(self, point: PointType):
        for i in range(1,len(self.listOuterVertex) - 1):
            if(NonConvexDomain.belongTriangle(point, [self.listOuterVertex[0],self.listOuterVertex[i],self.listOuterVertex[i+1]])):
                return True
        return False

    def addBoundaryPoint(self, point : PointType) -> int:
        if point not in self.listOuterVertex:
            for index, edge in enumerate(self.outerBoundary):
                if(NonConvexDomain.belongSegment(point, [self.listOuterVertex[edge[i]] for i in [0,1]])):
                    numPoint = len(self.listOuterVertex)
                    self.listOuterVertex.append(point)
                    self.outerBoundary = self.outerBoundary[:index] + [[edge[0],numPoint],[numPoint, edge[1]]] + self.outerBoundary[(index+1):]
                    return numPoint
            raise ValueError("The point given is not in the outer boundary.")
        else:
            for i, point2 in enumerate(self.listOuterVertex):
                if point2 == point:
                    return i

    def addWallPoint(self, point: PointType) -> int:
        if point not in self.wallVertices:
            for index, edge in enumerate(self.wallEdges):
                if(NonConvexDomain.belongSegment(point, [self.wallVertices[edge[i]] for i in [0,1]])):
                    numPoint = len(self.wallVertices)
                    self.wallVertices.append(point)
                    self.wallEdges = self.wallEdges[:index] + [[edge[0],numPoint],[numPoint, edge[1]]] + self.wallEdges[(index+1):]
                    return numPoint
            raise ValueError("The point given is not in an inner wall.")
        else:
            for i, point2 in enumerate(self.wallVertices):
                if point2 == point:
                    return i

    def addWall(self, coordWall: List[PointType], cycle: bool=False ) -> None:
        """
        Can only exclude convex domains.
        """
        for P in coordWall:
            if(P not in self):
                raise ValueError("The wall is not inside the domain.")

        if(cycle): #CircularWall
            for i in range(len(self.wallVertices),len(self.wallVertices)+len(coordWall)-1):
                self.wallEdges.append([i,i+1])
            self.wallEdges.append([len(self.wallVertices)+len(coordWall)-1, len(self.wallVertices)])

            for P in coordWall:
                self.wallVertices.append(P)

            Centerpoint = [0,0]
            for P in coordWall[:-1]:
                Centerpoint[0] += P[0]
                Centerpoint[1] += P[1]

            self.wallHolesPoint.append([Centerpoint[0]/(len(coordWall)-1),Centerpoint[1]/(len(coordWall)-1)])
        else:
            for i in range(len(self.wallVertices),len(self.wallVertices)+len(coordWall)-1):
                self.wallEdges.append([i,i+1])

            for P in coordWall:
                self.wallVertices.append(P)

    def addExit(self, exitEdge: List[PointType]):
        self.exitList.append(exitEdge)
        VertexIndex = [0,0]
        if(len(exitEdge) > 2):
            raise ValueError("An exit is a list of two points.")
        if(self.hasEdgeInInnerWall(exitEdge)):
            for i in [0,1]:
                VertexIndex[i] = self.addWallPoint(exitEdge[i])
            if VertexIndex not in self.wallEdges and [VertexIndex[1],VertexIndex[0]] not in self.wallEdges:
                raise ValueError("Impossible exit for this domain inner walls")
        elif(self.hasEdgeInOuterBoundary(exitEdge)):
            for i in [0,1]:
                VertexIndex[i] = self.addBoundaryPoint(exitEdge[i])
            if VertexIndex not in self.outerBoundary and [VertexIndex[1],VertexIndex[0]] not in self.outerBoundary:
                raise ValueError("Impossible exit for this domain boundary")
        else:
            raise ValueError("The exit is not inside an inner wall nor in an outer oboundary.")

    def addExits(self, ListOfExits: List[List[PointType]]):
        for exitEdge in ListOfExits:
            self.addExit(exitEdge)

    def findBoundaryPoint(self, P: PointType) -> int :
        for i,edge in enumerate(self.outerBoundary):
            if NonConvexDomain.belongSegment(P, [self.listOuterVertex[edge[0]], self.listOuterVertex[edge[1]]]):
                return i
        return -1

    def findWallPoint(self, P: PointType) -> int :
        for i, edge in enumerate(self.wallEdges):
            if NonConvexDomain.belongSegment(P, [self.wallVertices[edge[0]], self.wallVertices[edge[1]]]):
                return i
        return -1

    def findExitPoint(self, P: PointType) -> int:
        for i, exit in enumerate(self.exitList):
            if NonConvexDomain.belongSegment(P, exit):
                return i
        return -1

    def hasBoundaryPoint(self, P: PointType) -> bool :
        return (self.findBoundaryPoint(P) != -1)

    def hasWallPoint(self, P: PointType) -> bool :
        return (self.findWallPoint(P) != -1)

    def hasExitPoint(self, P: PointType) -> bool:
        return (self.findExitPoint(P) != -1)

    def hasExitEdge(self, edge : List[PointType]) -> bool:
        exitIndices = [self.findExitPoint(edge[0]), self.findExitPoint(edge[1])]
        if exitIndices[0] == -1 or exitIndices[1] == -1:
            return False
        if exitIndices[0] == exitIndices[1]:
            return True
        if NonConvexDomain.belongSegment(self.exitList[exitIndices[0]][0], self.exitList[exitIndices[1]]) or NonConvexDomain.belongSegment(self.exitList[exitIndices[0]][1], self.exitList[exitIndices[1]]):
            if np.abs( (self.exitList[exitIndices[0]][1][0]- self.exitList[exitIndices[0]][0][0])*(self.exitList[exitIndices[1]][1][1]- self.exitList[exitIndices[1]][0][1]) - (self.exitList[exitIndices[0]][1][1]- self.exitList[exitIndices[0]][0][1])*(self.exitList[exitIndices[1]][1][0]- self.exitList[exitIndices[1]][0][0]) ) < float(1e-10):
                return True
        return False

    def hasWallEdge(self, edge : List[PointType]) -> bool:
        wallIndices = [self.findWallPoint(edge[0]), self.findWallPoint(edge[1])]
        if wallIndices[0] == -1 or wallIndices[1] == -1:
            return False
        wall1Coord = [self.wallVertices[self.wallEdges[wallIndices[0]][0]],self.wallVertices[self.wallEdges[wallIndices[0]][1]]]
        wall2Coord = [self.wallVertices[self.wallEdges[wallIndices[1]][0]],self.wallVertices[self.wallEdges[wallIndices[1]][1]]]
        if NonConvexDomain.belongSegment(edge[0],wall2Coord) or NonConvexDomain.belongSegment(edge[1],wall1Coord):
            return True
        if NonConvexDomain.belongSegment(wall1Coord[0], wall2Coord) or NonConvexDomain.belongSegment(wall1Coord[1], wall2Coord):
            if np.abs( (wall1Coord[1][0]- wall1Coord[0][0])*(wall2Coord[1][1]- wall2Coord[0][1]) - (wall1Coord[1][1]- wall1Coord[0][1])*(wall2Coord[1][0]- wall2Coord[0][0]) ) < float(1e-10):
                return True
        return False

    def hasOuterEdge(self, edge : List[PointType]) -> bool:
        outerIndices = [self.findBoundaryPoint(edge[0]), self.findBoundaryPoint(edge[1])]
        if outerIndices[0] == -1 or outerIndices[1] == -1:
            return False
        outer1Coord = [self.listOuterVertex[self.outerBoundary[outerIndices[0]][0]],self.listOuterVertex[self.outerBoundary[outerIndices[0]][1]]]
        outer2Coord = [self.listOuterVertex[self.outerBoundary[outerIndices[1]][0]],self.listOuterVertex[self.outerBoundary[outerIndices[1]][1]]]
        if NonConvexDomain.belongSegment(edge[0],outer2Coord) or NonConvexDomain.belongSegment(edge[1],outer1Coord):
            return True
        if NonConvexDomain.belongSegment(outer1Coord[0], outer2Coord) or NonConvexDomain.belongSegment(outer1Coord[1], outer2Coord):
            if np.abs( (outer1Coord[1][0]- outer1Coord[0][0])*(outer2Coord[1][1]- outer2Coord[0][1]) - (outer1Coord[1][1]- outer1Coord[0][1])*(outer2Coord[1][0]- outer2Coord[0][0]) ) < float(1e-10):
                if(edge[1][0] == 4.6875):
                    print("outerIndices : ", outerIndices)
                    print("Oui 2 outers different mais alignées sur edge ", edge)
                    print("Resultat determinant : ",(outer1Coord[1][0]- outer1Coord[0][0])*(outer2Coord[1][1]- outer2Coord[0][1]) - (outer1Coord[1][1]- outer1Coord[0][1])*(outer2Coord[1][0]- outer2Coord[0][0]))
                    print("Avec Outer coords = ", outer1Coord, outer2Coord)
                return True
        return False

    def hasEdgeInOuterBoundary(self, edge: List[PointType]) -> bool:
        for i in [0,1]:
            if self.hasBoundaryPoint(edge[i]):
                return True
        return False

    def hasEdgeInInnerWall(self, edge: List[PointType]) -> bool:
        for i in [0,1]:
            if self.hasWallPoint(edge[i]):
                return True
        return False

    def show(self) -> None:
        if(plt):
            print("MatplotLib show not implemented yet.")
        elif(go):
            fig = go.Figure()
            self.addPlot(fig)
            fig.update_layout(yaxis=dict(
                scaleanchor='x',
                scaleratio=1))
            fig.show()
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")

    def addPlot(self, fig) -> None:
        if(go):
            startPoint = self.listOuterVertex[self.outerBoundary[0][0]]
            orderedOuterVertices = [startPoint]
            endPoint = self.listOuterVertex[self.outerBoundary[0][1]]
            while(endPoint != startPoint):
                orderedOuterVertices.append(endPoint)
                for edge in self.outerBoundary:
                    if(self.listOuterVertex[edge[0]] == endPoint):
                        endPoint = self.listOuterVertex[edge[1]]
                        break

            fig.add_trace(go.Scatter(x=[P[0] for P in orderedOuterVertices]+[orderedOuterVertices[0][0]],
                                    y=[P[1] for P in orderedOuterVertices]+[orderedOuterVertices[0][1]],
                                    fill="toself", fillcolor="White", mode="lines"))
            for edge in self.wallEdges:
                fig.add_shape(type="line",
                        x0=self.wallVertices[edge[0]][0],
                        y0=self.wallVertices[edge[0]][1],
                        x1=self.wallVertices[edge[1]][0],
                        y1=self.wallVertices[edge[1]][1],
                        line=dict(
                            color="LightSeaGreen",
                            width=2,
                        ))
            for exit in self.exitList:
                fig.add_shape(type="line",
                        x0=exit[0][0],
                        y0=exit[0][1],
                        x1=exit[1][0],
                        y1=exit[1][1],
                        line=dict(
                            color="Red",
                            width=2,
                        ))
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")

    @staticmethod
    def belongTriangle(M: PointType,T: List[PointType]) -> bool:
        det = (T[1][0]-T[0][0])*(T[2][1]-T[0][1]) - (T[1][1]-T[0][1])*(T[2][0]-T[0][0])
        if det == 0:
            return False

        X = ((M[0]-T[0][0])*(T[2][1]-T[0][1]) - (M[1]-T[0][1])*(T[2][0]-T[0][0]))/det
        Y = ((T[1][0]-T[0][0])*(M[1]-T[0][1]) - (T[1][1]-T[0][1])*(M[0]-T[0][0]))/det
        if( 0 <= X and 0 <= Y and X+Y <= 1):
            return True
        else:
            return False

    def belongSegment(P: PointType,AB: List[PointType]) -> bool:
        A = AB[0]
        B = AB[1]
        if(np.abs((A[0]- P[0])*(A[1]-B[1]) -(A[1]- P[1])*(A[0]-B[0])) < float(1e-10)):
            scal = (A[0] - P[0])*(A[0] - B[0]) + (A[1] - P[1])*(A[1] - B[1])
            if( scal <= (A[0] - B[0])**2 + (A[1] - B[1])**2 and scal >= 0):
                return True
        return False


class Mesh(object):
    """
    Mesh Object : object generating a triangular mesh of a given grain for a given domain object;
    It contains all the useful lists of edges, triangles and vertices (detailled in the init function)
    Attributes :
    - self.dx : float                                 # maximal area of a triangle in the mesh

    - self.vertices : ArrayLike                       # array of all the vertices TriangleEdgeCoordinates
    - self.edges : ArrayLike                          # array of all the edges as [vertexIndex, vertexIndex]
    - self.triangles : ArrayLike                      # array of all the triangles as [vertexIndex, vertexIndex, vertexIndex]

    - self.exitVertices : ArrayLike                   # array of the vertices index that belong to an exit.
    - self.exitEdges : ArrayLike                      # array of the edges index where the edge belongs to an exit
    - self.wallEdges : ArrayLike                      # array of the edges index where the edge belongs to an exit
    - self.boundaryPoints : list = []                 #
    - self.boundaryEdgesIndex : ArrayLike             #

    - self.trianglesWithEdges : ArrayLike             # array of triangles ordered as self.triangles, elements as [edgeIndex, edgeIndex, edgeIndex]
    - self.pairsOfTriangles : List[List[int]]         # nested list of length 1 or 2, ordered as self.edges, elements as [triangleIndex, triangleIndex] or [triangleIndex]
    - self.trianglesPerVertex : List[list]            # nested list, ordered as self.vertices, an element is a list of (number of triangles containing the vertex) elements as [[triangle index, [otherVertex1, otherVertex2 ]], ...]

    - self.outerNormalVectByTriangles : ArrayLike     #
    - self.cellAreas : ArrayLike                      #
    - self.edgeLength : ArrayLike                     #


    Methods :
    - get methods : classical get methods
    - compute methods : triangulate, computeEdgeList/EdgeLength/OuterNormals...
    correspond to all the redundant classical calculations related to the mesh
    - show :  plot the mesh
    """


    def __init__(self):
        self.dx : float = None

        self.vertices : ArrayLike
        self.edges : ArrayLike #Liste des edges par paires d'indices de vertex
        self.triangles : ArrayLike

        self.exitVertices : ArrayLike
        self.exitEdges : ArrayLike #Liste des edges d'exit par indice d'edge
        self.wallEdges : ArrayLike #Liste des wall edges par indice d'edge
        self.boundaryPoints : list = []
        self.boundaryEdgesIndex : ArrayLike #Liste des indices des edges faisant le bord du domaine
        self.vertexFlags : list = []

        self.trianglesWithEdges : ArrayLike #Liste des triangles par indice d'edge
        self.pairsOfTriangles : List[List[int]] #Liste des triangles entourant les edges ordonnée comme la EdgeList par les indices des triangles.
        self.trianglesPerVertex : List[list] #Nested List ordered by the vertex indices with the following content : [[triangle index, [otherVertex1, otherVertex2 ]], ...]

        self.outerNormalVectByTriangles : ArrayLike #Liste de triplets correspondants aux 3 vecteurs normaux unitaires dans le même ordre que les edges correspondant dans la TriangleWithEdgeList
        self.cellAreas : ArrayLike
        self.barycenters : ArrayLike
        self.edgeLength : ArrayLike

    def getEdgeCoord(self, edge):
        return [self.vertices[edge[0]], self.vertices[edge[1]]]

    def getNeighborPoints(self, index):
        listIndex = []
        for triangle in self.triangles:
            if(index in triangle):
                for j in triangle:
                    if(j != index and j not in listIndex):
                        listIndex.append(j)
        return(listIndex)

    def importMeshFromMsh(self, filename):
        if not meshio:
            raise ImportError("meshio must be installed in order to use .msh related methods.")
        Imesh = meshio.read(filename)
        self.vertices = np.array([ [P[0],P[1]] for P in Imesh.points ])
        self.triangles = Imesh.cells_dict["triangle"]

        self.computeEdgeList()
        self.exitEdges = []
        self.wallEdges = []

        print(Imesh)
        print(Imesh.cells_dict["line"])
        print(Imesh.point_data)
        #print(Imesh.point_data["gmsh:dim_tags"])
        print(Imesh.cell_data["gmsh:geometrical"])
        #print(Imesh.cell_data["gmsh:physical"])

        specialEdgesIndex = []
        if("line" in Imesh.cells_dict.keys()):
            for index, edge in enumerate(self.edges):
                for outerEdge in Imesh.cells_dict["line"]:
                    if (edge[0] == outerEdge[0] and edge[1] == outerEdge[1]) or (edge[0] == outerEdge[1] and edge[1] == outerEdge[0]) :
                        specialEdgesIndex.append(index)
                        break
            if(len(specialEdgesIndex) != len(Imesh.cells_dict["line"])):
                raise ImportError("Error importing .msh : some line cells are not edges from the mesh.")

        exitEdges = []
        wallEdges = []
        if("hughes2d:special" in Imesh.cell_data.keys()):
            for meshIndex, edgeIndex in enumerate(specialEdgesIndex):
                if Imesh.cell_data["hughes2d"][0][meshIndex] == 1:
                    exitEdges.append(edgeIndex)
                elif Imesh.cell_data["hughes2d"][0][meshIndex] == 2:
                    wallEdges.append(edgeIndex)
                else:
                    print("Warning : some special edges are neither wall nor exit.")
            self.exitEdges = np.array(exitEdges)
            self.wallEdges = np.array(wallEdges)

        else :
            print("WARNING : Mesh without walls imported from msh. Setting all special edges as walls.")
            for edgeIndex in specialEdgesIndex:
                if edgeIndex not in exitEdges:
                    wallEdges.append(edgeIndex)
            self.wallEdges = np.array(wallEdges)

        if(len(self.exitEdges) == 0):
            print("WARNING : Mesh without exits imported from msh.")
        if(len(specialEdgesIndex) == 0 and len(self.exitEdges) == 0 and len(self.wallEdges) == 0):
            raise ImportError("Error importing .msh : Mesh without walls neither exits")
        if(len(specialEdgesIndex) != len(self.exitEdges) + len(self.wallEdges)):
            print("WARNING : Losing some special edges during the import of .msh file.")

        print("Mesh imported. Contains %d triangles."% len(self.triangles))

        self.computePairOfTrianglesList()

        self.computeOuterNormals()

        self.dx = self.computeCellAreas()

        self.computeEdgeLength()

        self.computeTrianglesPerVertex()

    def importMeshFromMshFreeFem(self, filename : str, flag_dict : dict = {"domain" : 0, "exit" : 98, "wall" : 99}) -> None:
        with open(filename, "r") as file:
            line = file.readline().split()
            nb_vertices, nb_triangles, nb_spe_edges = int(line[0]), int(line[1]), int(line[2])
            print(nb_vertices, nb_triangles, nb_spe_edges)

            vertices = []
            exitVertices = []
            wallVertices = []
            for i in range(nb_vertices):
                line = file.readline().split()
                vertices.append([float(line[0]), float(line[1])])
                if(int(line[2]) == flag_dict['exit']):
                    exitVertices.append(i)
                elif(int(line[2]) == flag_dict['wall']):
                    wallVertices.append(i)
            self.vertices = np.array(vertices)
            self.exitVertices = np.array(exitVertices)
            self.boundaryPoints = np.array(exitVertices+wallVertices)

            triangles = []
            for i in range(nb_triangles):
                line = file.readline().split()
                triangles.append([int(line[0])-1, int(line[1])-1, int(line[2])-1])
            self.triangles = np.array(triangles)

            self.computeEdgeList()

            #print(self.triangles)

            exitEdges = []
            wallEdges = []

            for i in range(nb_spe_edges):
                line = file.readline().split()
                for index, edge in enumerate(self.edges):
                    if (edge[0] == int(line[0])-1 and edge[1] == int(line[1])-1 ) or (edge[0] == int(line[1])-1  and edge[1] == int(line[0])-1 ) :
                        if(int(line[2]) == flag_dict['exit']):
                            exitEdges.append(index)
                        elif(int(line[2]) == flag_dict['wall']):
                            wallEdges.append(index)
                        else :
                            print("Warning : some special edges are neither wall nor exit")
                        break

            self.exitEdges = np.array(exitEdges)
            self.wallEdges = np.array(wallEdges)


            if(len(self.exitEdges) == 0):
                print("WARNING : Mesh without exits imported from msh.")
            if(len(self.exitEdges) == 0 and len(self.wallEdges) == 0):
                raise ImportError("Error importing .msh : Mesh without walls neither exits")
            if(nb_spe_edges != len(self.exitEdges) + len(self.wallEdges)):
                print("WARNING : Losing some special edges during the import of .msh file.")

            print("Mesh imported. Contains %d triangles."% len(self.triangles))

            print(len(self.vertices), len(self.triangles), len(self.exitEdges)+len(self.wallEdges))

            print(len(self.edges), len(self.trianglesWithEdges))

            self.computePairOfTrianglesList()

            self.computeOuterNormals()

            self.dx = self.computeCellAreas()
            print("Maximal area for a triangle in the mesh : ", self.dx)

            self.computeEdgeLength()

            self.computeTrianglesPerVertex()

    def addExit(self, exitCoordinates):
        print("Pas imp")

    def addExits(self, exitList):
        print("pas imp")

    def addInnerWall(self, wall):
        print("pas imp")

    def addEdgeConstraint(self, edges, constraintLvl):
        print("pas imp")

    def exportMeshMsh(self, filename):
        if not meshio:
            raise ImportError("meshio must be installed in order to use .msh related methods.")
        points = self.vertices
        cells = [
            ("line", np.array([self.edges[i] for i in np.concatenate([self.exitEdges, self.wallEdges])])),
            ("triangle", self.triangles)
        ]



        cell_data_dict = {
        "gmsh:geometrical" : [[4 for i in np.concatenate([self.exitEdges, self.wallEdges])], np.array([1 for t in self.triangles])]
        #"gmsh:physical" : [[3 for i in np.concatenate([self.exitEdges, self.wallEdges])], np.array([1 for t in self.triangles])],
        #"hughes2d:special" : [np.concatenate([[1 for i in self.exitEdges],[2 for w in self.wallEdges]]), np.array([0 for t in self.triangles])]
        #"exits" : [self.edges[i] for i in self.exitEdges],
        #"walls" : [self.edges[i] for i in self.wallEdges]
        }

        mesh = meshio.Mesh(
            points,
            cells,
            point_data = {"gmsh:dim_tags" : np.array([[2,0] for P in self.vertices])},
            # Optionally provide extra data on points, cells, etc.
            cell_data = cell_data_dict
            )
        mesh.write(
            filename,  # str, os.PathLike, or buffer/open file
            file_format="gmsh"  # optional if first argument is a path; inferred from extension
        )
        print("Saving mesh as ", filename)

    def exportMeshMshFreeFem(self, filename):
        with open(filename, "w") as file:
            file.write("%d %d %d\n"% (len(self.vertices), len(self.triangles), len(self.exitEdges)+len(self.wallEdges)))
            #file.write("Vertices")
            for i,point in enumerate(self.vertices):
                file.write("%f %f %d\n"% (point[0], point[1], self.vertexFlags[i]))
            #file.write("Triangles")
            for triangle in self.triangles:
                file.write("%d %d %d 0\n"%(triangle[0]+1,triangle[1]+1,triangle[2]+1))
            #file.write("Edges")
            for edgeIndex in self.exitEdges:
                file.write("%d %d 98\n"%(self.edges[edgeIndex][0]+1,self.edges[edgeIndex][1]+1))
            for edgeIndex in self.wallEdges:
                file.write("%d %d 99\n"%(self.edges[edgeIndex][0]+1,self.edges[edgeIndex][1]+1))

        print("Saving mesh as ", filename)

    def getDomainFromMesh(self) -> NonConvexDomain:
        print("pas imp")
        return None

    def generateMeshFromDomain(self, domain: NonConvexDomain, dx: float, da: float=30) -> None:
        """
        Compute a triangular mesh covering domain with the maximal length of an edge being set to dx.
        The computations are done using triangle
        da : min angle insiade a triangle (be careful, crash if set too high)
        dx : max area of a triangle (heavy computations when set low, computation time ~ 1/(dx)^2 )
        """
        domainVertices = []
        domainSpecialEdges = []
        domainVertices += domain.listOuterVertex
        domainSpecialEdges += domain.outerBoundary

        N = len(domainVertices)
        domainVertices += domain.wallVertices
        for edge in domain.wallEdges:
            domainSpecialEdges.append([edge[0]+N, edge[1] +N])

        if(len(domain.wallHolesPoint)):
            D = dict(vertices=domainVertices, segments=domainSpecialEdges,holes=domain.wallHolesPoint)
        else:
            D = dict(vertices=domainVertices, segments=domainSpecialEdges)

        flag = 'q'+str(da)+'pa'+str(dx)
        meshDict = tr.triangulate(D, flag)

        if 'segments' not in meshDict.keys():
            raise ValueError("Corrupted domain, impossible to generate a mesh.")
        print("Mesh generated. Contains %d triangles" % len(triangles))

        boundaryPoints = []
        for edge in meshDict['segments']:
            for index in edge:
                if(index not in boundaryPoints):
                    boundaryPoints.append(index)

        self.boundaryPoints = np.array(boundaryPoints)

        self.vertices = vertices
        self.triangles = triangles

        self.computeEdgeList()

        self.setExitsFromDomain(domain)

        self.computeVertexFlags(domain)

        #print(self.wallEdges)
        #print(np.array([[self.vertices[self.edges[edge][0]],self.vertices[self.edges[edge][1]]] for edge in self.wallEdges]))
        #BUG à fix hasExitEdge Domain dans les angles....

        self.computePairOfTrianglesList()

        self.computeOuterNormals()

        self.dx = self.computeCellAreas()

        self.computeEdgeLength()

        self.computeTrianglesPerVertex()


    def computeEdgeList(self) -> None:
        """
        Fills two arrays :
        - self.edges : cointaining all the edges of the mesh as a pair of vertex index.
        - self.trianglesWithEdges : containing all the triangles of self.triangles, in the same order, represented as a triplet of edge index.

        These arrays are required for : ......
        """
        edges = [] #Liste des edges par paires d'indices de vertex
        TriangleWithEdgeList = [] #Liste des triangles par indice d'edge

        for triangle in self.triangles:
            Edge1 = [triangle[0],triangle[1]]
            Edge2 = [triangle[1],triangle[2]]
            Edge3 = [triangle[2],triangle[0]]

            L = [Edge1,Edge2,Edge3]
            TriangleEdgeCoordinates = []

            for edge2 in L:
                isAlreadyIn = 0
                for n,edge in enumerate(edges):
                    if( edge[0] == edge2[0] and edge[1] == edge2[1] ) or ( edge[0] == edge2[1] and edge[1] == edge2[0] ):
                        isAlreadyIn = 1
                        TriangleEdgeCoordinates.append(n)
                        break
                if( not isAlreadyIn):
                    edges.append(edge2)
                    TriangleEdgeCoordinates.append(len(edges)-1)

            TriangleWithEdgeList.append(TriangleEdgeCoordinates)
        self.edges = np.array(edges)
        self.trianglesWithEdges = np.array(TriangleWithEdgeList)

    def setExitsFromDomain(self, domain: NonConvexDomain) -> None:
        """
        Configure the exits and the wall edges and vertices lists from the domain object passed in parameter.

        INPUT :
        - domain : NonConvexDomain instance

        These lists are required for : .....
        """
        exitVertices = []
        for i in range(len(self.vertices)):
            for exit in domain.exitList:
                if NonConvexDomain.belongSegment(self.vertices[i],exit):
                    exitVertices.append(i)
        self.exitVertices = np.array(exitVertices)
        #print("Exit Vertices : ", self.ExitVertices)

        exitEdges = [] # List of Edge exit par indices d'edge
        for i,edge in enumerate(self.edges):
            if(domain.hasExitEdge([self.vertices[edge[0]],self.vertices[edge[1]]])):
                exitEdges.append(i)
        self.exitEdges = np.array(exitEdges)

        if(len(exitEdges) == 0):
            raise ValueError("Your mesh has no exit edge.")

        wallEdges = []
        #Juste remplir avec les methodes de non convex domain....
        for index, edge in enumerate(self.edges):
                if( ( domain.hasWallEdge([self.vertices[edge[0]],self.vertices[edge[1]]]) and not domain.hasExitEdge([self.vertices[edge[0]],self.vertices[edge[1]]]) )
                        or ( domain.hasOuterEdge([self.vertices[edge[0]],self.vertices[edge[1]]]) and not domain.hasExitEdge([self.vertices[edge[0]],self.vertices[edge[1]]]) )):
                    wallEdges.append(index)
        self.wallEdges = np.array(wallEdges)

    def computeVertexFlags(self, domain : NonConvexDomain, flag_dict : dict = {"domain" : 0, "exit" : 98, "wall" : 99}) -> None:
        self.vertexFlags = []
        for point in self.vertices:
            if domain.hasExitPoint(point):
                self.vertexFlags.append(flag_dict["exit"])
            elif domain.hasWallPoint(point) or domain.hasBoundaryPoint(point):
                self.vertexFlags.append(flag_dict["wall"])
            else:
                self.vertexFlags.append(flag_dict["domain"])


    def computeTrianglesPerVertex(self) -> None:
        """
        Compute the self.trianglesPerVertex list.

        Required for : EikonalSolver, ....
        """
        self.trianglesPerVertex : list = []

        for index in range(len(self.vertices)):
            L :list = []
            for triangleIndex, triangle in enumerate(self.triangles):
                if index in triangle:
                    Exclude : list = [int(P) for P in triangle if P != index]
                    L.append([triangleIndex, Exclude])
            self.trianglesPerVertex.append(L)


    def computePairOfTrianglesList(self) -> None:
        """
        Compute the lists pairsOfTriangles and boundaryEdgesIndex;

        Required for : .......
        """
        self.pairsOfTriangles=[] #Liste des triangles entourant les edges ordonnée comme la EdgeList par les indices des triangles.
        boundaryEdgesIndex = [] #Liste des indices des edges faisant le bord du domaine

        for index, edge in enumerate(self.edges):
            PairOfTriangles = []
            NumberOfLoop = 0
            while(len(PairOfTriangles) < 2 and NumberOfLoop < len(self.trianglesWithEdges)):
                if(index in self.trianglesWithEdges[NumberOfLoop]):
                    PairOfTriangles.append(NumberOfLoop)
                NumberOfLoop += 1
            self.pairsOfTriangles.append(PairOfTriangles)
            if(len(PairOfTriangles) < 2):
                boundaryEdgesIndex.append(index)
        self.boundaryEdgesIndex = np.array(boundaryEdgesIndex)

    def computeOuterNormals(self) -> None:
        outerNormalVectByTriangles = []
        for index, triangle in enumerate(self.trianglesWithEdges):
            L = []
            for edgeindex in triangle:
                for Pindex in self.triangles[index]:
                    if(Pindex not in self.edges[edgeindex]):
                        L.append(ComputeOuterNormalUnitVect(self.vertices[self.edges[edgeindex][0]], self.vertices[self.edges[edgeindex][1]], self.vertices[Pindex]))
            if(len(L) != 3):
                print("erreur")
                return
            outerNormalVectByTriangles.append(L)
        self.outerNormalVectByTriangles = np.array(outerNormalVectByTriangles)

    def computeCellAreas(self) -> float:
        """
        Compute the areas of each triangular cell. Raises an error if a triangle is degenerated.
        Computes also the barycenter of each triangle.

        Return :
        - the maximal area of a triangle in the mesh

        """
        barycenters = []
        cellAreas = []
        max = 0
        for triangle in self.triangles:
            A = self.vertices[triangle[0]]
            B = self.vertices[triangle[1]]
            C = self.vertices[triangle[2]]

            barycenters.append([(A[0]+B[0]+C[0])/3,(A[1]+B[1]+C[1])/3])
            cellAreas.append(abs( (B[0]-A[0])*(C[1]-A[1]) - (B[1]-A[1])*(C[0]-A[0]))/2)
            if(cellAreas[-1] > max):
                max = cellAreas[-1]
            elif(cellAreas[-1] == 0):
                print(A,B,C)
                raise ValueError("Degenerated mesh : at least one of the triangle is of area 0.")
        self.cellAreas = np.array(cellAreas)
        self.barycenters = np.array(barycenters)
        return(max)

    def computeEdgeLength(self) -> None:
        """
        Compute the edgeLength array.
        """
        edgeLength = []
        for edge in self.edges:
            A = self.vertices[edge[0]]
            B = self.vertices[edge[1]]

            edgeLength.append(np.sqrt( (B[0]-A[0])*(B[0]-A[0]) + (B[1]-A[1])*(B[1]-A[1])))
        self.edgeLength = np.array(edgeLength)

    def appendDict(self, dico):

        MeshDico = {"type": "triangular mesh"}
        MeshDico["dx"] = self.dx
        MeshDico["vertices"] = self.vertices.tolist()
        MeshDico["edges"] = self.edges
        MeshDico["triangles"] = self.triangles.tolist()

        MeshDico["listOuterVertex"] = [self.vertices[i].tolist() for i in self.boundaryPoints]

        MeshDico["exitVertices"] = self.exitVertices.tolist()
        MeshDico["exitEdges"] = self.exitEdges.tolist() #Liste des edges d'exit par indice d'edge
        MeshDico["wallEdges"] = self.wallEdges.tolist() #Liste des wall edges par indice d'edge
        MeshDico["boundaryPoints"] = self.boundaryPoints.tolist()
        MeshDico["vertexFlags"] = self.vertexFlags.tolist() #Liste des triangles entourant les edges ordonnée comme la EdgeList par les indices des triangles.
        MeshDico["boundaryEdgesIndex"] = self.boundaryEdgesIndex.tolist() #Liste des indices des edges faisant le bord du domaine

        MeshDico["trianglesWithEdges"] = self.trianglesWithEdges.tolist()
        MeshDico["pairsOfTriangles"] = self.pairsOfTriangles
        MeshDico["trianglesPerVertex"] = self.trianglesPerVertex

        MeshDico["outerNormalVectByTriangles"] = self.outerNormalVectByTriangles.tolist() #Liste de triplets correspondants aux 3 vecteurs normaux unitaires dans le même ordre que les edges correspondant dans la TriangleWithEdgeList
        MeshDico["cellAreas"] = self.cellAreas.tolist()
        MeshDico["barycenters"] = self.barycenters.tolist()
        MeshDico["edgeLength"] = self.edgeLength.tolist()

        dico["mesh"] = MeshDico

    def saveToJson(self, filename):
        if not json:
            raise ImportError("Module json not or wrongly installed. Needed for the json methods")

        MeshDico = {"type": "triangular mesh"}
        MeshDico["dx"] = self.dx
        MeshDico["vertices"] = self.vertices.tolist()
        MeshDico["edges"] = self.edges.tolist()
        MeshDico["triangles"] = self.triangles.tolist()

        MeshDico["listOuterVertex"] = [self.vertices[i].tolist() for i in self.boundaryPoints]

        MeshDico["exitVertices"] = self.exitVertices.tolist()
        MeshDico["exitEdges"] = self.exitEdges.tolist() #Liste des edges d'exit par indice d'edge
        MeshDico["wallEdges"] = self.wallEdges.tolist() #Liste des wall edges par indice d'edge
        MeshDico["boundaryPoints"] = self.boundaryPoints.tolist()
        MeshDico["vertexFlags"] = self.vertexFlags #Liste des triangles entourant les edges ordonnée comme la EdgeList par les indices des triangles.
        MeshDico["boundaryEdgesIndex"] = self.boundaryEdgesIndex.tolist() #Liste des indices des edges faisant le bord du domaine

        MeshDico["trianglesWithEdges"] = self.trianglesWithEdges.tolist()
        MeshDico["pairsOfTriangles"] = self.pairsOfTriangles
        MeshDico["trianglesPerVertex"] = self.trianglesPerVertex

        MeshDico["outerNormalVectByTriangles"] = self.outerNormalVectByTriangles.tolist() #Liste de triplets correspondants aux 3 vecteurs normaux unitaires dans le même ordre que les edges correspondant dans la TriangleWithEdgeList
        MeshDico["cellAreas"] = self.cellAreas.tolist()
        MeshDico["barycenters"] = self.barycenters.tolist()
        MeshDico["edgeLength"] = self.edgeLength.tolist()

        with open(filename+"_mesh.json", 'w', encoding='utf-8') as f:
            json.dump(MeshDico, f, ensure_ascii=False, indent=4)

    def loadFromJson(self, filename):
        if not json:
            raise ImportError("Module json not or wrongly installed. Needed for the json methods")
        with open(filename) as f:
            data = json.load(f)
            self.dx : float = data['dx']

            self.vertices : ArrayLike = np.array(data['vertices'])
            self.edges : ArrayLike = np.array(data['edges']) #Liste des edges par paires d'indices de vertex
            self.triangles : ArrayLike = np.array(data['triangles'])

            self.exitVertices : ArrayLike = np.array(data['exitVertices'])
            self.exitEdges : ArrayLike = np.array(data['exitEdges'])#Liste des edges d'exit par indice d'edge
            self.wallEdges : ArrayLike = np.array(data['wallEdges'])#Liste des wall edges par indice d'edge
            self.boundaryPoints : list = data['boundaryPoints']
            self.boundaryEdgesIndex : ArrayLike = np.array(data['boundaryEdgesIndex'])#Liste des indices des edges faisant le bord du domaine
            self.vertexFlags : list = data['vertexFlags']

            self.trianglesWithEdges : ArrayLike = np.array(data['trianglesWithEdges'])#Liste des triangles par indice d'edge
            self.pairsOfTriangles : List[List[int]] = data['pairsOfTriangles'] #Liste des triangles entourant les edges ordonnée comme la EdgeList par les indices des triangles.
            self.trianglesPerVertex : List[list] = data['trianglesPerVertex'] #Nested List ordered by the vertex indices with the following content : [[triangle index, [otherVertex1, otherVertex2 ]], ...]

            self.outerNormalVectByTriangles : ArrayLike = np.array(data['outerNormalVectByTriangles'])#Liste de triplets correspondants aux 3 vecteurs normaux unitaires dans le même ordre que les edges correspondant dans la TriangleWithEdgeList
            self.cellAreas : ArrayLike = np.array(data['cellAreas'])
            self.barycenters : ArrayLike = np.array(data['barycenters'])
            self.edgeLength : ArrayLike = np.array(data['edgeLength'])

            print("Mesh successfully loaded.")



    def show(self, with_domain: NonConvexDomain=None) -> None:
        if(go):
            fig = go.Figure()
            if(with_domain):
                with_domain.addPlot(fig)
            for T in self.triangles:
                fig.add_trace(go.Scatter(x=[self.vertices[i][0] for i in T]+[self.vertices[T[0]][0]],
                                        y=[self.vertices[i][1] for i in T]+[self.vertices[T[0]][1]],
                                fill="toself",
                                fillcolor="White",
                                mode="lines",
                                line=dict(
                                    color="Black",
                                    width=1
                                 )))
            fig.update_layout(yaxis=dict(
                scaleanchor='x',
                scaleratio=1))
            for edge in self.wallEdges:
                fig.add_shape(type="line",
                        x0=self.vertices[self.edges[edge][0]][0],
                        y0=self.vertices[self.edges[edge][0]][1],
                        x1=self.vertices[self.edges[edge][1]][0],
                        y1=self.vertices[self.edges[edge][1]][1],
                        line=dict(
                            color="LightSeaGreen",
                            width=2,
                        ))
            for edge in self.exitEdges:
                fig.add_shape(type="line",
                        x0=self.vertices[self.edges[edge][0]][0],
                        y0=self.vertices[self.edges[edge][0]][1],
                        x1=self.vertices[self.edges[edge][1]][0],
                        y1=self.vertices[self.edges[edge][1]][1],
                        line=dict(
                            color="Red",
                            width=2,
                        ))

            fig.show()
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")

    def addPlot(self, fig):
        if(go):
            for T in self.triangles:
                fig.add_trace(go.Scatter(x=[self.vertices[i][0] for i in T]+[self.vertices[T[0]][0]],
                                        y=[self.vertices[i][1] for i in T]+[self.vertices[T[0]][1]],
                                fill="toself",
                                fillcolor="White",
                                mode="lines",
                                line=dict(
                                    color="Black",
                                    width=1
                                 )))
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")


class CellValueMap(object):
    """
    Cell Value map :
    correspond to a function which is constant on the triangles (ex : densities)
    """

    def __init__(self, Mesh):
        self.Mesh = Mesh
        self.values = [0 for _ in self.Mesh.triangles]

    def generateRandom(self,variability = 0.4):
        self.values = [ 0.2 + variability*(alea.random()-0.5) for _ in self.Mesh.triangles]

    def __len__(self):
        return len(self.values)

    def __getitem__(self, index): #renvoie le tuple (triangle, value)
        return self.values[index]

    def __setitem__(self, index, value): #renvoie le tuple (triangle, value)
        self.values[index] = value

    def __add__(self, other):
        result = CellValueMap(self.Mesh)
        for index in range(len(self.Mesh.triangles)):
            result[index] = self[index] + other[index]
        return result

    def __iadd__(self, other):
        for index in range(len(self.Mesh.triangles)):
            self[index] += other[index]

    def __mul__(self,a):
        result = CellValueMap(self.Mesh)
        for index in range(len(self.Mesh.triangles)):
            result[index] = self[index]*a
        return result

    def __imul__(self, a):
        for index in range(len(self.Mesh.triangles)):
            self[index] *= a

    def valueOnTriangle(self, L):
        for index, triangle in enumerate(self.Mesh.triangles):
            if(list(triangle) == list(L)):
                return self.values[index]
        return None

    def setConstantCircle(self, center:List[PointType], radius:float, value:float):
        for index in range(len(self.Mesh.triangles)):
            if((self.Mesh.barycenters[index][0] - center[0])**2 + (self.Mesh.barycenters[index][1] - center[1])**2 <= radius):
                self.values[index] = value
            else:
                self.values[index] = 0

    def fitAveragedMap(self, other): #Realise la moyenne pondérée entre deux CellValueMap dont le mesh est potentiellement différent.
        for i in range(len(self.values)):
            LenCell = self.Mesh.points[i+1] - self.Mesh.points[i]
            start = 0
            while(other.Mesh.points[start] < self.Mesh.points[i]):
                start += 1
            if(other.Mesh.points[start] > self.Mesh.points[i]):
                start -= 1
            end = start+1
            while(other.Mesh.points[end] < self.Mesh.points[i+1]):
                end += 1
            j = start
            Average = 0
            while(j <= end -1):
                Average += other.values[j]*(min(other.Mesh.points[j+1], self.Mesh.points[i+1]) - max(other.Mesh.points[j], self.Mesh.points[i]))
                j += 1
            self.values[i] = Average/LenCell

    def show(self):
        if(go):
            fig = go.Figure()
            #self.Mesh.domain.addPlot(fig)
            for j,T in enumerate(self.Mesh.triangles):
                fig.add_trace(go.Scatter(x=[self.Mesh.vertices[i][0] for i in T]+[self.Mesh.vertices[T[0]][0]],
                                        y=[self.Mesh.vertices[i][1] for i in T]+[self.Mesh.vertices[T[0]][1]],
                                fill="toself",
                                hoverinfo = "none",
                                showlegend = False,
                                mode="none",
                                fillcolor ='rgb('+str( int(255*min(1,max(self.values[j],0))) )+',0,0)'
                                ))
            fig.update_layout(yaxis=dict(
                scaleanchor='x',
                scaleratio=1))
            fig.show()
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")

    def getScatter(self):
        if(go):
            L = []
            for j,T in enumerate(self.Mesh.triangles):
                L.append(go.Scatter(x=[self.Mesh.vertices[i][0] for i in T]+[self.Mesh.vertices[T[0]][0]],
                                        y=[self.Mesh.vertices[i][1] for i in T]+[self.Mesh.vertices[T[0]][1]],
                                fill="toself",
                                hoverinfo = "none",
                                showlegend = False,
                                mode="none",
                                fillcolor ='rgb('+str( int(255*min(1,max(self.values[j],0))) )+',0,0)'
                                ))
            return L
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")





class VertexValueMap(object):
    """
    VertexValueMap
    corresponds to a function which is affine on the triangles, defined by its values on the vertices
    Methods :
    - computeGradientFlow : return a list of normalized gradients by triangle.
    - show : show the vertices values
    - show vector fields : show the normalized gradients
    """

    def __init__(self, Mesh):
        self.Mesh = Mesh
        self.values = [0 for _ in self.Mesh.vertices]

    def generateRandom(self,variability = 0.5):
        self.values = [ 0.5 + variability*(alea.random()-0.5) for _ in self.Mesh.vertices]

    def __len__(self):
        return len(self.values)

    def __getitem__(self, index): #renvoie le tuple (triangle, value)
        return self.values[index]

    def __setitem__(self,index, value):
        self.values[index] = value

    def __str__(self):
        return(str(self.values))

    def setInfinity(self):
        self.values = [float('inf') for _ in self.Mesh.vertices]

    def checkGradientValidity(self):
        Grads = self.computeGradientFlow()
        for edgeIndex in self.Mesh.WallEdges:
            indexTriangle = self.Mesh.pairsOfTriangles[edgeIndex][0]
            Grad = Grads[indexTriangle]
            for P in self.Mesh.triangles[indexTriangle]:
                if P not in self.Mesh.EdgeList[edgeIndex]:
                    OuterVect = ComputeOuterNormalUnitVect(self.Mesh.vertices[self.Mesh.EdgeList[edgeIndex][0]],
                                                            self.Mesh.vertices[self.Mesh.EdgeList[edgeIndex][1]],
                                                            self.Mesh.vertices[P])
                    C = P
            Scal = OuterVect[0]*Grad[0] + OuterVect[1]*Grad[1]
            if(Scal > 0):
                print("Vecteur non valide ! Emplacement du edge : ")
                print(self.Mesh.vertices[self.Mesh.EdgeList[edgeIndex][0]], self.Mesh.vertices[self.Mesh.EdgeList[edgeIndex][1]])
                print("Valeurs : ")
                print(self.values[self.Mesh.EdgeList[edgeIndex][0]],self.values[self.Mesh.EdgeList[edgeIndex][1]], self.values[C])
                print("Intensité de l'erreur : ")
                print(Scal)

    def show3D(self):
        if(go):
            fig = go.Figure()
            fig.add_trace(go.Mesh3d(x=[P[0] for P in self.Mesh.vertices],
                                    y = [P[1] for P in self.Mesh.vertices],
                                    z = self.values,
                                    opacity=1,
                                    color='rgba(244,22,100,0.6)'))
            fig.show()
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")

    def show(self, grid=False, colorscale_name = 'viridis'):
        if(go):
            fig = go.Figure()
            #self.Mesh.domain.addPlot(fig)
            if(grid):
                self.Mesh.addPlot(fig)
            fig.add_trace(go.Scatter(x=[P[0] for P in self.Mesh.vertices],
                                    y = [P[1] for P in self.Mesh.vertices],
                            hoverinfo = "none",
                            showlegend = False,
                            mode="markers",
                            marker = dict(
                            color = self.values,
                            colorscale = colorscale_name
                            )))
            fig.update_layout(yaxis=dict(
                scaleanchor='x',
                scaleratio=1))
            fig.show()
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")


    def computeGradientFlow(self):
        if(go):
            LTrianglesGrad = []
            for point in range(len(self.Mesh.vertices)):
                if(self.values[point] < 0):
                    print("GROS PROBLEME");
            for triangle in self.Mesh.triangles: #On a des déterminants égaux à 0....
                    det =  ((self.Mesh.vertices[triangle[1]][0] - self.Mesh.vertices[triangle[0]][0])*(self.Mesh.vertices[triangle[2]][1] - self.Mesh.vertices[triangle[0]][1])
                            - (self.Mesh.vertices[triangle[2]][0] - self.Mesh.vertices[triangle[0]][0])*(self.Mesh.vertices[triangle[1]][1] - self.Mesh.vertices[triangle[0]][1]))
                    if(det == 0):
                        print("det nul :")
                        print(str(self.Mesh.vertices[triangle[0]][0])+','+str(self.Mesh.vertices[triangle[0]][1]))
                        print(str(self.Mesh.vertices[triangle[1]][0])+','+str(self.Mesh.vertices[triangle[1]][1]))
                        print(str(self.Mesh.vertices[triangle[2]][0])+','+str(self.Mesh.vertices[triangle[2]][1]))
                    Vecx = ( (self.values[triangle[0]] - self.values[triangle[2]])*(self.Mesh.vertices[triangle[1]][1] - self.Mesh.vertices[triangle[0]][1])
                            + (self.values[triangle[1]] - self.values[triangle[0]])*(self.Mesh.vertices[triangle[2]][1] - self.Mesh.vertices[triangle[0]][1]) )/det
                    Vecy = ( (self.values[triangle[0]] - self.values[triangle[1]])*(self.Mesh.vertices[triangle[2]][0] - self.Mesh.vertices[triangle[0]][0])
                            + (self.values[triangle[2]] - self.values[triangle[0]])*(self.Mesh.vertices[triangle[1]][0] - self.Mesh.vertices[triangle[0]][0]) )/det
                    LTrianglesGrad.append([-Vecx/np.sqrt(Vecx**2 + Vecy**2),-Vecy/np.sqrt(Vecx**2 + Vecy**2)])
            return LTrianglesGrad
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")

    def showVectorField(self):
        if(go):
            L = self.computeGradientFlow()
            fig = go.Figure()
            #self.Mesh.domain.addPlot(fig)
            fig.add_trace(go.Scatter(x=[P[0] for P in self.Mesh.vertices],
                                    y = [P[1] for P in self.Mesh.vertices],
                            hoverinfo = "none",
                            showlegend = False,
                            mode="markers",
                            marker = dict(
                            color = self.values,
                            colorscale = 'viridis'
                            )))
            figQuiv = ff.create_quiver([(self.Mesh.vertices[T[0]][0]+self.Mesh.vertices[T[1]][0]+self.Mesh.vertices[T[2]][0])/3 for T in self.Mesh.triangles],
                                        [(self.Mesh.vertices[T[0]][1]+self.Mesh.vertices[T[1]][1]+self.Mesh.vertices[T[2]][1])/3 for T in self.Mesh.triangles],
                                        [V[0] for V in L], [V[1] for V in L])
            fig.add_traces(figQuiv.data)
            fig.update_layout(yaxis=dict(
                scaleanchor='x',
                scaleratio=1))
            fig.show()
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")




def ComputeOuterNormalUnitVect(A,B,C):
    N = np.sqrt((A[0]-B[0])**2 + (A[1]-B[1])**2)
    x = (N**2 + (C[0]-B[0])**2 + (C[1]-B[1])**2 - (C[0]-A[0])**2 - (C[1]-A[1])**2 )/(2*N) #Distance AH
    """if(x < 0):
        print("triangle non aigu...") #Point d'intersection entre AB et la hauteur issue de C"""
    H = [ A[0] + (B[0]-A[0])*x/N, A[1] + (B[1]-A[1])*x/N ]
    Sign = (C[0]-H[0])*(B[1]-A[1]) + (C[1]-H[1])*(A[0]-B[0])
    if(Sign > 0):
        return [(A[1]-B[1])/N,(B[0]-A[0])/N]
    return [(B[1]-A[1])/N,(A[0]-B[0])/N]
