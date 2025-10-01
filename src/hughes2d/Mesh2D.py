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
    import matplotlib.patches as patches
    import matplotlib.cm as cm
    import matplotlib.collections as collections
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
    """Object corresponding to the geometry of a domain. The domain can have walls and exits on its boundary as well as walls inside the domain.

    Examples:
        Creating a square domain::

            MyDomain = NonConvexDomain([[0,0],[0,3],[3,3],[3,0]])

        Adding an exit::

            MyDomain.addExit([[0,1],[0,2]])

        Creating a hole in the center of the domain::

            MyDomain.addWall([[1,1],[1,2],[2,2],[2,1]], cycle=True)

    Args:
        outerVerticesList (List[List[float]): The ordered list of the coordinates of the vertices defining the boundary of the domain. A vertex is represented as a pair of floats i.e. [float,float].

    Attributes
    -----------

    Attributes:
        outerVertices (List[List[float]]): The list of the coordinates of the vertices located on the boundary of the domain. A vertex is represented as a pair of floats i.e. [float,float].
        wallVertices (List[List[float]]): The list of the coordinates of the vertices located on walls inside the domain. A vertex is represented as a pair of floats i.e. [float,float].
        outerBoundary (List[List[int]]): The list of the edges of the outer boundary of the domain. An edge is represented as a pair of indices of the vertices of outerVertices i.e. [int,int].
        wallEdges (List[List[int]]): The list of the edges of the outer boundary of the domain. An edge is represented as a pair of indices of the vertices of outerVertices i.e. [int,int].
        wallHolesPoint (List[List[float]]): The list of the coordinates of points located inside the holes of the domains. A point is represented as a pair of floats i.e. [float,float].
        exitList (List[int]): The list of the edges corresponding to exits of the domain. An edge is represented as a pair of indices of the vertices of outerVertices i.e. [int,int].
        zones (dict):

    Methods
    --------
    """

    def __init__(self, outerVerticesList: List[PointType]=[[0,0],[1,0],[0,1]] ) -> None:
        self.outerVertices = outerVerticesList
        self.outerBoundary = []
        for i in range(len(outerVerticesList)):
            self.outerBoundary.append([i,(i+1)%(len(outerVerticesList))])
        #print(self.outerBoundary)
        self.wallVertices = []
        self.wallEdges = []
        self.wallHolesPoint = []
        self.exitList = []
        self.zones = dict()

    def importFromDXF(self, filename: str) -> int:
        """
        Import from a .dxf file. The dxf file can contain three layers :

        - a "domain" layer (mandatory) that contain the outer boundary
        - an "innerWalls" layer that contains all the inner structures of the domain.
        - an "exits" layer that contains all the exits of the domain.
        - (optional) "zone_<zone_name>" layers defining the zones with zone name.

        The package ezdxf must be installed in order to use this method.

        Args:
            filename (str): the path to the .dxf file to import.

        Returns:
            int: Error code
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

            self.outerVertices = []
            self.outerBoundary = []
            self.wallVertices = []
            self.wallEdges = []
            self.wallHolesPoint = []

            # iterate over all entities in modelspace
            msp = doc.modelspace()

            # entity query for all LINE entities in modelspace
            for e in msp.query("LINE"):
                if(e.dxf.layer == "0" or e.dxf.layer == "domain"):
                    P = [np.round(e.dxf.start[0],5),np.round(e.dxf.start[1],5)]
                    if P not in self.outerVertices:
                        numP = len(self.outerVertices)
                        self.outerVertices.append(P)
                    else:
                        numP = self.addBoundaryPoint(P)

                    Q = [np.round(e.dxf.end[0],5),np.round(e.dxf.end[1],5)]
                    if Q not in self.outerVertices:
                        numQ = len(self.outerVertices)
                        self.outerVertices.append(Q)
                    else:
                        numQ = self.addBoundaryPoint(Q)

                    self.outerBoundary.append([numP, numQ])

            if([len(self.outerVertices)-1, 0] not in self.outerBoundary ):
                raise ValueError("Corrupted domain : domain not closed.")
            wallVert = []
            for e in msp.query("LINE"):
                if(e.dxf.layer == "innerWalls"):
                    cycling = 0
                    P = [np.round(e.dxf.start[0],5),np.round(e.dxf.start[1],5)]
                    if P not in wallVert:
                        wallVert.append(P)
                    else:
                        cycling += 1

                    Q = [np.round(e.dxf.end[0],5),np.round(e.dxf.end[1],5)]
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
                        print("Adding wall : ", wallVert)
                        self.addWall(wallVert, cycle =True)
                        wallVert = []

            if(len(wallVert) > 0):
                self.addWall(wallVert, cycle =False)

            for e in msp.query("LINE"):
                if(e.dxf.layer == "exits"):
                    P = [np.round(e.dxf.start[0],5),np.round(e.dxf.start[1],5)]
                    Q = [np.round(e.dxf.end[0],5),np.round(e.dxf.end[1],5)]
                    self.addExit([P,Q])

            zoneDict = dict()
            for e in msp.query("LINE"):
                if(str(e.dxf.layer).find("zone_") == 0):
                    label = str(e.dxf.layer)[5:]
                    if(label not in zoneDict.keys()):
                        zoneDict[label] = []

                    P = [np.round(e.dxf.start[0],5),np.round(e.dxf.start[1],5)]
                    if P not in zoneDict[label]:
                        zoneDict[label].append(P)

                    Q = [np.round(e.dxf.end[0],5),np.round(e.dxf.end[1],5)]
                    if Q not in zoneDict[label]:
                        zoneDict[label].append(Q)

            for label in zoneDict.keys():
                self.addZone(label, zoneDict[label])

            print(self.zones)
        else:
            raise ImportError("ezdxf not or wrongly installed.")


    def __contains__(self, point: PointType):
        """
        This method tests if the point passed as a parameter is inside the convex hull of the domain.

        Args:
            point (List[float]): The coordinates of the point to test respresented as [float, float].

        Return:
            bool
        """
        for i in range(1,len(self.outerVertices) - 1):
            if(belongTriangle(point, [self.outerVertices[0],self.outerVertices[i],self.outerVertices[i+1]])):
                return True
        return False

    def addBoundaryPoint(self, point : PointType) -> int:
        """
        Add a given point of the boundary to the list of points of the domain.
        This is typically used in order to guarantee that the given point will be included as a vertex of the mesh generated from this domain.
        If the point is already in outerVertices, the method returns the index of the point without adding it.

        Args:
            point (List[float]): The point to add to the boundary respresented as [float, float].

        Raises:
            ValueError: The method raises a value error if the given point is not in the boundary of the domain.

        Returns:
            int: The method returns the index corresponding to the added point in the outerVertices list.
        """
        if point not in self.outerVertices:
            for index, edge in enumerate(self.outerBoundary):
                if(NonConvexDomain.belongSegment(point, [self.outerVertices[edge[i]] for i in [0,1]])):
                    numPoint = len(self.outerVertices)
                    self.outerVertices.append(point)
                    self.outerBoundary = self.outerBoundary[:index] + [[edge[0],numPoint],[numPoint, edge[1]]] + self.outerBoundary[(index+1):]
                    return numPoint
            raise ValueError("The point given is not in the outer boundary.")
        else:
            for i, point2 in enumerate(self.outerVertices):
                if point2 == point:
                    return i

    def getLimits(self):
        """
        Computes and return the extremal bounds of the domain

        Returns:
            List[List[float]]: the extremal coordinates of the domain as a list of two two-elements lists i.e. [[x_min,x_max],[y_min,y_max]].
        """
        x_min, x_max = self.outerVertices[0][0],self.outerVertices[0][0]
        y_min, y_max = self.outerVertices[0][1],self.outerVertices[0][1]
        for P in self.outerVertices:
            if P[0] > x_max:
                x_max = P[0]
            if P[1] > y_max:
                y_max = P[1]
            if P[0] < x_min:
                x_min = P[0]
            if P[1] < y_min:
                y_min = P[1]
        return [[x_min,x_max],[y_min,y_max]]


    def addZone(self, zoneName:str, zoneVertices:List[PointType]) -> None:
        """
        Adds a zone to the dict of zones of the domain.

        Args:
            zoneName (str): the name to be given to the new zone.
            zoneVertices (List[List[float]]): list of the vertex of the boundary of the zone.

        Raises:
            NameError: raises a name error if the zone name is already taken.
        """
        if(zoneName in self.zones.keys()):
            raise NameError("The zone name is already used.")
        self.zones[zoneName] = zoneVertices


    def addWallPoint(self, point: PointType) -> int:
        """
        Add a given point of the boundary to the list of points of the inner walls of the domain.
        This is typically used in order to guarantee that the given point will be included as a vertex of the mesh generated from this domain.
        If the point is already in wallVertices, the method returns the index of the point without adding it.

        Args:
            point (List[float]): the point to add on an inner wall of the domain.

        Raises:
            ValueError: the method raises a value error if the given point is not in the inner walls of the domain.

        Returns:
            int: the method returns the index corresponding to the added point in the wallVertices list.
        """
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
        Adds a wall in the domain. If cycle is set to True, the wall is considered as an area to exclude from the domain.
        If not, only the edges defined by the coordWall are excluded from the domain.

        Note:
            The cycle=True works properly only with convex holes in the domain. In fact, the barycenters of all the points of the walls to add must be inside the hole for the exclusion to work properly.

        Args:
            coordWall (List[PointType]) : List of the vertices defining the walls.
            cycle (bool) : A boolean switching between a hole to exclude from the domain and only walls of zero thickness to exclude from the domain.

        Raises:
            ValueError: raises a Value Error if the points of wallCoord are not inside the convex hull of the domain.
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
            for P in coordWall:
                Centerpoint[0] += P[0]
                Centerpoint[1] += P[1]

            self.wallHolesPoint.append([Centerpoint[0]/len(coordWall),Centerpoint[1]/len(coordWall)])
        else:
            for i in range(len(self.wallVertices),len(self.wallVertices)+len(coordWall)-1):
                self.wallEdges.append([i,i+1])

            for P in coordWall:
                self.wallVertices.append(P)

    def addExit(self, exitEdge: List[PointType]):
        """
        Adds an exit passed in parameter. If the two exit points are on different edges, the shortest path (in number of vertices crossed) following the wall is added as an exit.

        Args:
        - exitEdge (List[PointType]) : an exit defined by two extremal points.

        Raises:
            ValueError: raises an error if:

                - at least one of the extremal points is not in a wall edge or a boundary edge;
                - there exists no path staying either in the walls or in the boundary that links the two extremal points of the exit.

        """
        VertexIndex = [0,0]
        if(len(exitEdge) > 2):
            raise ValueError("An exit is a list of two points.")
        if(self.hasWallEdge(exitEdge)):
            for i in [0,1]:
                VertexIndex[i] = self.addWallPoint(exitEdge[i])

            exitPath = NonConvexDomain.shortestPathBFS(VertexIndex[0],VertexIndex[1],self.wallEdges)
            if len(exitPath) == 0:
                raise ValueError("Impossible exit for this domain inner walls")
            else:
                exitEdges = [ [self.wallVertices[exitPath[i]],self.wallVertices[exitPath[i+1]]] for i in range(len(exitPath) - 1)]
                self.exitList += exitEdges
        elif(self.hasOuterEdge(exitEdge)):
            for i in [0,1]:
                VertexIndex[i] = self.addBoundaryPoint(exitEdge[i])
            exitPath = NonConvexDomain.shortestPathBFS(VertexIndex[0],VertexIndex[1],self.outerBoundary)
            if len(exitPath) == 0:
                raise ValueError("Impossible exit for this domain boundary")
            else:
                exitEdges = [ [self.outerVertices[exitPath[i]],self.outerVertices[exitPath[i+1]]] for i in range(len(exitPath) - 1)]
                self.exitList += exitEdges
        else:
            raise ValueError("The exit is not inside an inner wall nor in an outer oboundary.")

    def addExits(self, ListOfExits: List[List[PointType]]):
        """
        Call the method addExit multiple times.

        Args:
            ListOfExits (List[List[PointType]]): a list of the exits to add.
        """
        for exitEdge in ListOfExits:
            self.addExit(exitEdge)

    def findBoundaryPoint(self, point: PointType) -> int :
        """
        Searches an edge of the boundary containing P.

        Args:
            point (PointType): the point the should be in the researched edge.

        Returns:
            int: the index of the first edge found in the boundary that contains P; returns -1 if no such edge is found.
        """
        for i,edge in enumerate(self.outerBoundary):
            if NonConvexDomain.belongSegment(point, [self.outerVertices[edge[0]], self.outerVertices[edge[1]]]):
                return i
        return -1

    def findWallPoint(self, point: PointType) -> int :
        """
        Searches an edge of the inner walls containing a given point.

        Args:
            point (List[float]): the point the should be in the researched edge.

        Returns:
            int: the index of the first edge found in the inner walls that contains the point; returns -1 if no such edge is found.
        """
        for i, edge in enumerate(self.wallEdges):
            if NonConvexDomain.belongSegment(point, [self.wallVertices[edge[0]], self.wallVertices[edge[1]]]):
                return i
        return -1

    def findExitPoint(self, point: PointType) -> int:
        """
        Searches an edge of the exits containing a given point.

        Args:
            point (List[float]): the point the should be in the researched edge.

        Returns:
            int: the index of the first edge found in the exits that contains the point; returns -1 if no such edge is found.
        """
        for i, exit in enumerate(self.exitList):
            if NonConvexDomain.belongSegment(point, exit):
                return i
        return -1

    def hasBoundaryPoint(self, point: PointType) -> bool :
        """
        Checks that the point passed as a parameter belongs to the outer boundary.
        Args:
            point (List[float]): the point to test.

        Returns:
            bool: A boolean that is True if the point is in the outer boundary.
        """
        return (self.findBoundaryPoint(point) != -1)

    def hasWallPoint(self, point: PointType) -> bool :
        """
        Checks that the point passed as a parameter belongs to the inner walls.
        Args:
            point (List[float]): the point to test.

        Returns:
            bool: A boolean that is True if the point is in the inner walls.
        """
        return (self.findWallPoint(point) != -1)

    def hasExitPoint(self, point: PointType) -> bool:
        """
        Checks that the point passed as a parameter belongs to the exits.
        Args:
            point (List[float]): the point to test.

        Returns:
            bool: A boolean that is True if the point is in the exits.
        """
        return (self.findExitPoint(point) != -1)

    def hasExitEdge(self, edge : List[PointType]) -> bool:
        """
        Tests if the segment defined by the pair of points passed as a parameter is inside an exit.

        Args:
            edge (List[PointType]) :  a pair of points defining a segment to test.

        Returns:
            bool
        """
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
        """
        Tests if the segment defined by the pair of points passed as a parameter is inside an inner wall.

        Args:
            edge (List[PointType]) :  a pair of points defining a segment to test.

        Returns:
            bool
        """
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

        verticesIndex = [0,0]
        for index, point in enumerate(self.wallVertices):
            for i in range(2):
                if ((edge[i][0] - point[0])**2 + (edge[i][1] - point[1])**2) < float(1e-10):
                    verticesIndex[i] = index
        for wall in self.wallEdges:
            if( wall[0] == verticesIndex[0] and wall[1] == verticesIndex[1]) or ( wall[0] == verticesIndex[1] and wall[1] == verticesIndex[0]):
                return True

        return False

    def hasOuterEdge(self, edge : List[PointType]) -> bool:
        """
        Tests if the segment defined by the pair of points passed as a parameter is inside the outer boundary.

        Args:
            edge (List[PointType]) :  a pair of points defining a segment to test.

        Returns:
            bool
        """
        outerIndices = [self.findBoundaryPoint(edge[0]), self.findBoundaryPoint(edge[1])]
        if outerIndices[0] == -1 or outerIndices[1] == -1:
            return False
        outer1Coord = [self.outerVertices[self.outerBoundary[outerIndices[0]][0]],self.outerVertices[self.outerBoundary[outerIndices[0]][1]]]
        outer2Coord = [self.outerVertices[self.outerBoundary[outerIndices[1]][0]],self.outerVertices[self.outerBoundary[outerIndices[1]][1]]]
        if NonConvexDomain.belongSegment(edge[0],outer2Coord) or NonConvexDomain.belongSegment(edge[1],outer1Coord):
            return True
        if NonConvexDomain.belongSegment(outer1Coord[0], outer2Coord) or NonConvexDomain.belongSegment(outer1Coord[1], outer2Coord):
            if np.abs( (outer1Coord[1][0]- outer1Coord[0][0])*(outer2Coord[1][1]- outer2Coord[0][1]) - (outer1Coord[1][1]- outer1Coord[0][1])*(outer2Coord[1][0]- outer2Coord[0][0]) ) < float(1e-10):
                return True
        return False

    def show(self, preference ="plotly") -> None:
        """
        Plotting method for the domain object. The method opens either a new window or a default navigator tab.

        Args:
            preference (str, optional): set to "plotly" or "matplotlib" to chose the preferred plotting package. If only one package is installed the preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed for this method to work. If both are installed and no preference is set, plotly is used.
        """
        if( go and (not plt or preference == "plotly") ):
            fig = go.Figure()
            self.addPlot(fig)
            fig.update_layout(yaxis=dict(
                scaleanchor='x',
                scaleratio=1))
            fig.show()
        elif(plt):
            fig, ax = plt.subplots()
            self.addPlot(fig,ax,preference)
            limits = self.getLimits()
            ax.set_xlim(limits[0][0],limits[0][1])
            ax.set_ylim(limits[1][0],limits[1][1])
            plt.axis('equal')
            plt.show()
            fig.clf()
            plt.close(fig)
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")

    def addPlot(self, fig, ax=None, preference="plotly") -> None:
        """
        Non-blocking plotting method for the domain object. The method does not show the graph.

        Args:
            fig (matplotlib.Figure or plotly.graph_objects.Figure): the instance of Figure to which the plot must be added. Depends on the library used.
            ax (matplotlib.Axes or None): the instance of Axes to use if matplotlib is used.
            preference (str): set to "plotly" or "matplotlib" to chose the preferred plotting package. If only one package is installed the preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed for this method to work. If both are installed and no preference is set, plotly is used.
        """
        if( go and (not plt or preference == "plotly") ): #plotly version of the plot
            startPoint = self.outerVertices[self.outerBoundary[0][0]]
            orderedOuterVertices = [startPoint]
            endPoint = self.outerVertices[self.outerBoundary[0][1]]
            while(endPoint != startPoint):
                orderedOuterVertices.append(endPoint)
                for edge in self.outerBoundary:
                    if(self.outerVertices[edge[0]] == endPoint):
                        endPoint = self.outerVertices[edge[1]]
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
        elif(plt): #matplotlib version of the plot
            startPoint = self.outerVertices[self.outerBoundary[0][0]]
            orderedOuterVertices = [startPoint]
            endPoint = self.outerVertices[self.outerBoundary[0][1]]
            while(endPoint != startPoint):
                orderedOuterVertices.append(endPoint)
                for edge in self.outerBoundary:
                    if(self.outerVertices[edge[0]] == endPoint):
                        endPoint = self.outerVertices[edge[1]]
                        break

            domain_polygon = patches.Polygon([(P[0],P[1]) for P in orderedOuterVertices], edgecolor='black', facecolor='white')
            ax.add_patch(domain_polygon)

            for edge in self.wallEdges:
                path = [    (self.wallVertices[edge[0]][0], self.wallVertices[edge[0]][1]),
                            (self.wallVertices[edge[1]][0], self.wallVertices[edge[1]][1]) ]
                ax.add_patch(patches.Polygon(path, edgecolor='black',linewidth=1))

            for exit in self.exitList:
                path = [    (exit[0][0], exit[0][1]),
                            (exit[1][0], exit[1][1]) ]
                ax.add_patch(patches.Polygon(path, edgecolor='red',linewidth=2))

        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")

    @staticmethod
    def shortestPathBFS(Istart: int , Iend: int, listEdges: List[List[int]]):
        """
        Computes the shortest path between Istart and Iend (in number of vertices) on the network defined by the integers as vertice and the edges of listEdges.
        The method uses a Breadth First Search algorithm with memory of the pathes explored.

        Args:
            Istart (int) : the starting vertex of the searched path.
            Iend (int) : the ending vertex of the searched path.
            listEdges (List[List[int]]) : a list of the edges of the network symbolized by a list of two integers.

        Returns:
            List[int]: The shortest path found as a list of integers i.e. the list of the successive vertices defining the path. If no path is found the method returns an empty list.
        """
        visited = []
        Pathes = [[Istart]]
        while(len(Pathes) > 0):
            newPathes = []
            for path in Pathes:
                if path[-1] == Iend:
                    return path
                visited.append(path[-1])

                for edge in listEdges:
                    if path[-1] == edge[0] and edge[1] not in visited:
                        newPathes.append(path+[edge[1]])
                    if path[-1] == edge[1] and edge[0] not in visited:
                        newPathes.append(path+[edge[0]])
            Pathes = newPathes

        return []

    def belongSegment(P: PointType,AB: List[PointType]) -> bool:
        """
        Checks that a point belongs to a segment within an error margin of float(1e-10).

        Args:
            P (PointType) : the point to test.
            AB (List[PointType]) : the segment to test.

        Returns:
            bool
        """
        A = AB[0]
        B = AB[1]
        if(np.abs((A[0]- P[0])*(A[1]-B[1]) -(A[1]- P[1])*(A[0]-B[0])) < float(1e-10)):
            scal = (A[0] - P[0])*(A[0] - B[0]) + (A[1] - P[1])*(A[1] - B[1])
            if( scal <= (A[0] - B[0])**2 + (A[1] - B[1])**2 and scal >= 0):
                return True
        return False

def belongTriangle(M: PointType,T: List[PointType]) -> bool:
    """
    Checks that a point belongs to a triangle.

    Args:
        M (PointType) : the point to test.
        T (List[PointType]) : the triangle to test.

    Returns:
        bool
    """
    det = (T[1][0]-T[0][0])*(T[2][1]-T[0][1]) - (T[1][1]-T[0][1])*(T[2][0]-T[0][0])
    if det == 0:
        return False

    X = ((M[0]-T[0][0])*(T[2][1]-T[0][1]) - (M[1]-T[0][1])*(T[2][0]-T[0][0]))/det
    Y = ((T[1][0]-T[0][0])*(M[1]-T[0][1]) - (T[1][1]-T[0][1])*(M[0]-T[0][0]))/det
    if( 0 <= X and 0 <= Y and X+Y <= 1):
        return True
    else:
        return False


class Mesh(object):
    """
    The Mesh Object is an object representing a triangular mesh of a given grain for a given domain object;
    It contains all the useful lists of edges, triangles and vertices for optimized computations of the numerical scheme involved in the hughes2d package.

    Examples:
        The mesh can be generated by using a NonConvexDomain as the basis of computations::

            MyDomain = NonConvexDomain([[0,0],[0,1],[1,1],[1,0]])
            MyMesh = Mesh()
            MyMesh.generateMeshFromDomain(MyDomain,0.1)

        It is also possible to import a mesh from a .msh file::

            MyMesh = Mesh()
            MyMesh.importMeshFromMsh("filename.msh")

        See also *importMeshFromMshFreeFem* and *importFromLists* for different import methods. The mesh can then be exported (see *exportMeshMsh* and *exportMeshMshFreeFem*) or saved as a .json file::

            MyMesh.saveToJson("filename")

        And can be loaded from this .json file::

            MyMesh2 = Mesh()
            MyMesh2.loadFromJson("filename.json")

    Attributes
    -----------

    Attributes :
        dx (float): maximal area of a triangle in the mesh
        vertices (ArrayLike): array of all the vertices TriangleEdgeCoordinates
        edges (ArrayLike): array of all the edges as [vertexIndex, vertexIndex]
        triangles (ArrayLike): array of all the triangles as [vertexIndex, vertexIndex, vertexIndex]
        exitVertices (ArrayLike): array of the vertices index that belong to an exit.
        exitEdges (ArrayLike): array of the edges index where the edge belongs to an exit
        wallEdges (ArrayLike): array of the edges index where the edge belongs to a wall
        boundaryPoints (ArrayLike): array of the indices of the vertices that belong to the outer boundary of the mesh
        boundaryEdgesIndex (ArrayLike): array of the edges index where the edge belongs to the outer boundary
        trianglesWithEdges (ArrayLike): array of triangles ordered as self.triangles, elements as [edgeIndex, edgeIndex, edgeIndex]
        pairsOfTriangles (List[List[int]]): nested list of length 1 or 2, ordered as self.edges, elements as [triangleIndex, triangleIndex] or [triangleIndex]
        trianglesPerVertex (List[list]): nested list, ordered as self.vertices, an element is a list of (number of triangles containing the vertex) elements as [[triangle index, [otherVertex1, otherVertex2 ]], ...]
        outerNormalVectByTriangles (ArrayLike): array order as self.tringles containing the unit normal vectors corresponding to the three edges of each triangle. The normal vectors are directed towards the exterior of  the triangle.
        cellAreas (ArrayLike): array ordered as self.triangles containing the area of each triangle of the mesh.
        edgeLength (ArrayLike):array ordered as self.edges containing the length of each edge.
        zones (dict):


    Methods
    --------
    """


    def __init__(self):
        self.dx : float = None

        self.vertices : ArrayLike
        self.edges : ArrayLike #Liste des edges par paires d'indices de vertex
        self.triangles : ArrayLike

        self.exitVertices : ArrayLike
        self.exitEdges : ArrayLike #Liste des edges d'exit par indice d'edge
        self.wallEdges : ArrayLike #Liste des wall edges par indice d'edge
        self.boundaryPoints : ArrayLike
        self.boundaryEdgesIndex : ArrayLike #Liste des indices des edges faisant le bord du domaine
        self.vertexFlags : list = []

        self.trianglesWithEdges : ArrayLike #Liste des triangles par indice d'edge
        self.pairsOfTriangles : List[List[int]] #Liste des triangles entourant les edges ordonnée comme la EdgeList par les indices des triangles.
        self.trianglesPerVertex : List[list] #Nested List ordered by the vertex indices with the following content : [[triangle index, [otherVertex1, otherVertex2 ]], ...]

        self.outerNormalVectByTriangles : ArrayLike #Liste de triplets correspondants aux 3 vecteurs normaux unitaires dans le même ordre que les edges correspondant dans la TriangleWithEdgeList
        self.cellAreas : ArrayLike
        self.barycenters : ArrayLike
        self.edgeLength : ArrayLike

        self.zones : dict = dict()

    def importMeshFromMsh(self, filename :str,verbose : bool = False, requirements : List[str] = ['all']):
        """
        Imports the data from a .msh file into the Mesh object.

        Args:
            filename (string) : the path to the file to import.
            verbose (bool, optional).
            requirements (List[str], optional): a list containing the computations that will be done using the mesh. The possible values are:

                - EikonalSolver : the mesh will be used in order to solve an eikonal equation
                - LWRSolver : the mesh will be used in order to solve a scalar conservation law
                - all : all the possible computations will be done
                - integrate : the mesh will be used in order to compute some integrals over the domain
                - FreeFEM : the mesh will be used in order to import/export mesh file for FreeFEM

        Raises:
            ImportError: requires the python library meshio.
        """
        if not meshio:
            raise ImportError("meshio must be installed in order to use .msh related methods.")
        Imesh = meshio.read(filename)
        self.vertices = np.array([ [P[0],P[1]] for P in Imesh.points ])
        self.triangles = Imesh.cells_dict["triangle"]

        self.computeEdgeList()
        self.exitEdges = []
        self.wallEdges = []

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
        self.computations(verbose=verbose, requirements=requirements)

    def importMeshFromMshFreeFem(self, filename : str, flag_dict : dict = {"domain" : 0, "exit" : 98, "wall" : 99}, verbose : bool = False, requirements : List[str] = ['all']) -> None:
        """
        Imports the data from a mesh file constructed in FreeFEM into the Mesh object.
        The specific structure (inner walls, exits...) can be specified by specifying different flags in FreeFEM.
        Args:
        - filename (str) : the path to the file to import.
        - flag_dict (dict) : a dictionary describing the specific translation of the FreeFEM flag number.
                    Must contain the keys domain, exit and wall.
        - verbose (bool, optional).
        - requirements (List[str], optional): a list containing the computations that will be done using the mesh. The possible values are:
                - EikonalSolver : the mesh will be used in order to solve an eikonal equation
                - LWRSolver : the mesh will be used in order to solve a scalar conservation law
                - all : all the possible computations will be done
                - integrate : the mesh will be used in order to compute some integrals over the domain
                - FreeFEM : the mesh will be used in order to import/export mesh file for FreeFEM
        """
        with open(filename, "r") as file:
            line = file.readline().split()
            nb_vertices, nb_triangles, nb_spe_edges = int(line[0]), int(line[1]), int(line[2])

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

            if verbose :
                print("Mesh imported. Contains %d triangles."% len(self.triangles))

            self.computations(verbose=verbose, requirements=requirements)

    def importFromLists(self, vertices:List[PointType], triangles:List[List[int]], domain:NonConvexDomain, verbose :bool = False, requirements : List[str] = ['all']):
        """
        Imports the lists passed as parameters into the Mesh object.

        Args:
            vertices (List[PointType]) : list of all the vertices coordinates.
            triangles (List[List[int]]) : list of the triangles symbolized by a triplet of the indices of the corresponding verttices in the vertices list.
            domain (NonConvexDomain) : a NonConvexDomain instance corresponding to the mesh imported.
            verbose (bool, optional): displaying progression in console.
            requirements (List[str], optional): a list containing the computations that will be done using the mesh. The possible values are:

                - EikonalSolver : the mesh will be used in order to solve an eikonal equation
                - LWRSolver : the mesh will be used in order to solve a scalar conservation law
                - all : all the possible computations will be done
                - integrate : the mesh will be used in order to compute some integrals over the domain
                - FreeFEM : the mesh will be used in order to import/export mesh file for FreeFEM
        """
        self.vertices = np.array(vertices)
        self.triangles = np.array(triangles)
        self.computations(domain, verbose, requirements)

    def addConvexZone(self, zoneName: str, zoneVertices: List[PointType]):
        """
        Adds a convex zone to the Mesh object.

        Args:
            zoneName (str) : the unique name designing the zone to add.
            zoneVertices (List[PointType]) :  the vertices of the outer boundary of the zone. The zone is supposed convex.
        """
        if(zoneName in self.zones.keys()):
            raise NameError("Zone name already in use.")

        self.zones[zoneName] = {'boundary' : zoneVertices,
                                'triangles' : []}

    def computeZonesTriangles(self):
        """
        Computes the triangles included in each zone.
        """
        for index, center in enumerate(self.barycenters):
            for zoneName in self.zones.keys():
                if self.inZone(center, zoneName):
                    self.zones[zoneName]['triangles'].append(index)

    def inZone(self, point, zoneName):
        """
        Tests if a given point is included in a given zone.

        Args:
            point (PointType) :  the point to test.
            zoneName (str) :  the name of the zone to test.
        """
        if zoneName not in self.zones.keys():
            raise NameError("The name "+zoneName+" does not correspond to a zone of the Mesh.")
        for i in range(1,len(self.zones[zoneName]['boundary']) - 1):
            if(belongTriangle(point, [self.zones[zoneName]['boundary'][0],self.zones[zoneName]['boundary'][i],self.zones[zoneName]['boundary'][i+1]])):
                return True
        return False

    def exportMeshMsh(self, filename):
        """
        Exports the data of the Mesh object in a .msh file.

        Args:
            filename (str) :  the complete path to the file.

        Raises:
            ImportError: requires meshio to be installed to work properly.
        """
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
        """
        Exports the data of the Mesh object in a .msh file following the specific structure of the FreeFEM mesh files.

        Args:
            filename (str) :  the complete path to the file.
        """
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


    def generateMeshFromDomain(self, domain: NonConvexDomain, dx: float, da: float=30, verbose :bool = False, requirements : List[str] = ['all']) -> None:
        """
        Compute a triangular mesh covering domain with the maximal length of an edge being set to dx.
        The computations are done using triangle (see https://www.cs.cmu.edu/~quake/triangle.html)

        Args:
            domain (NonConvexDomain): the domain object from which the mesh must be generated.
            dx (float): max area of a triangle (heavy computations when set low, computation time ~ 1/(dx)^2 ).
            da (float, optional): min angle inside a triangle (be careful, crash if set too high). Default = 30.
            verbose (bool, optional): displaying progression in console.
            requirements (List[str], optional): a list containing the computations that will be done using the mesh. The possible values are:

                - EikonalSolver : the mesh will be used in order to solve an eikonal equation
                - LWRSolver : the mesh will be used in order to solve a scalar conservation law
                - all : all the possible computations will be done
                - integrate : the mesh will be used in order to compute some integrals over the domain
                - FreeFEM : the mesh will be used in order to import/export mesh file for FreeFEM
        """
        domainVertices = []
        domainSpecialEdges = []
        domainVertices += domain.outerVertices
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
        if verbose :
            print("Mesh generated. Contains %d triangles" % len(meshDict['triangles']))

        boundaryPoints = []
        for edge in meshDict['segments']:
            for index in edge:
                if(index not in boundaryPoints):
                    boundaryPoints.append(index)

        self.boundaryPoints = np.array(boundaryPoints)

        self.vertices = meshDict['vertices']
        self.triangles = meshDict['triangles']
        self.computations(domain, verbose, requirements)

    def computations(self, domain: NonConvexDomain = None, verbose : bool = False, requirements : List[str] = ['all']) -> None:
        """
        Triggers all the computations required for the Mesh object.
        Args:
            domain (NonConvexDomain, optional): the domain object to use in order to recover the exits and the zones.
            verbose (bool, optional): displaying progression in console.
            requirements (List[str], optional): a list containing the computations that will be done using the mesh. The possible values are:

                - EikonalSolver : the mesh will be used in order to solve an eikonal equation
                - LWRSolver : the mesh will be used in order to solve a scalar conservation law
                - all : all the possible computations will be done
                - integrate : the mesh will be used in order to compute some integrals over the domain
                - FreeFEM : the mesh will be used in order to import/export mesh file for FreeFEM
        """
        if "all" in requirements:
            self.computeEdgeList()

            if domain is not None:
                self.setExitsFromDomain(domain)
                self.computeVertexFlags(domain)
                for zoneName in domain.zones.keys():
                    self.addConvexZone(zoneName, domain.zones[zoneName])

            self.computePairOfTrianglesList()
            self.computeOuterNormals()
            self.dx = self.computeCellAreas()
            if(verbose):
                print("Minimal area for a triangle in the mesh : ", self.dx)
            self.computeEdgeLength()
            self.computeTrianglesPerVertex()
            self.computeZonesTriangles()
        else:
            if("EikonalSolver" in requirements or "LWRSolver" in requirements or "FreeFEM" in requirements):
                self.computeEdgeList()

            if("EikonalSolver" in requirements or "LWRSolver" in requirements):
                self.setExitsFromDomain(domain)
                for zoneName in domain.zones.keys():
                    self.addConvexZone(zoneName, domain.zones[zoneName])
                self.computeEdgeLength()

            if("FreeFEM" in requirements):
                self.computeVertexFlags(domain)

            if("LWRSolver" in requirements):
                self.computePairOfTrianglesList()
                self.computeOuterNormals()

            if("LWRSolver" in requirements or "integrate" in requirements):
                self.dx = self.computeCellAreas()
                if(verbose):
                    print("Minimal area for a triangle in the mesh : ", self.dx)

            if("EikonalSolver" in requirements or "integrate" in requirements):
                self.computeTrianglesPerVertex()

            if("integrate" in requirements):
                self.computeZonesTriangles()


    def computeEdgeList(self) -> None:
        """
        Fills two arrays :

        - self.edges : cointaining all the edges of the mesh as a pair of vertex index.
        - self.trianglesWithEdges : containing all the triangles of self.triangles, in the same order, represented as a triplet of edge index.

        Note:
            Required for:

                - LWRSolver
                - EikonalSolver
                - setExitsFromDomain
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

        Args:
            domain (NonConvexDomain): the domain object containing information about the exits.

        Note:
            Required for:

                - LWRSolver
                - EikonalSolver
        """
        exitVertices = []
        for i in range(len(self.vertices)):
            for exit in domain.exitList:
                if NonConvexDomain.belongSegment(self.vertices[i],exit):
                    exitVertices.append(i)
        self.exitVertices = np.array(exitVertices)

        exitEdges = []
        for i,edge in enumerate(self.edges):
            if(domain.hasExitEdge([self.vertices[edge[0]],self.vertices[edge[1]]])):
                exitEdges.append(i)
        self.exitEdges = np.array(exitEdges)

        if(len(exitEdges) == 0):
            raise ValueError("Your mesh has no exit edge.")

        wallEdges = []
        for index, edge in enumerate(self.edges):
                if( ( domain.hasWallEdge([self.vertices[edge[0]],self.vertices[edge[1]]]) and not domain.hasExitEdge([self.vertices[edge[0]],self.vertices[edge[1]]]) )
                        or ( domain.hasOuterEdge([self.vertices[edge[0]],self.vertices[edge[1]]]) and not domain.hasExitEdge([self.vertices[edge[0]],self.vertices[edge[1]]]) )):
                    wallEdges.append(index)
        self.wallEdges = np.array(wallEdges)

    def computeVertexFlags(self, domain : NonConvexDomain, flag_dict : dict = {"domain" : 0, "exit" : 98, "wall" : 99}) -> None:
        """
        Computation of the vertices flags necessary to export the Mesh as FreeFEM mesh file.

        Args:
            domain (NonConvexDomain, optional): the domain object to use in order to recover the exits and the zones.
            flag_dict (dict): dictionary prescribing the integer corresponding to the different types of vertex for an export as a FreeFEM .msh file.

        Note:
            Required for:

                - exportMeshMshFreeFem
        """
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

        Note:
            Required for:

                - EikonalSolver
                - integrateOverSquareBall
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
        Compute the lists pairsOfTriangles and boundaryEdgesIndex.

        Note:
            Required for:

                - LWRSolver
                - checkGradientValidity
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
        """
        Compute the list outerNormalVectByTriangles.

        Note:
            Required for:

                - LWRSolver
        """

        outerNormalVectByTriangles = []

        def ComputeOuterNormalUnitVect(A,B,C):
            N = np.sqrt((A[0]-B[0])**2 + (A[1]-B[1])**2)
            x = (N**2 + (C[0]-B[0])**2 + (C[1]-B[1])**2 - (C[0]-A[0])**2 - (C[1]-A[1])**2 )/(2*N)
            H = [ A[0] + (B[0]-A[0])*x/N, A[1] + (B[1]-A[1])*x/N ]
            Sign = (C[0]-H[0])*(B[1]-A[1]) + (C[1]-H[1])*(A[0]-B[0])
            if(Sign > 0):
                return [(A[1]-B[1])/N,(B[0]-A[0])/N]
            return [(B[1]-A[1])/N,(A[0]-B[0])/N]
            
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

        Returns:
            float: the minimal area of a triangle in the mesh.

        Note:
            Required for:

                - integrate
                - LWRSolver
        """
        barycenters = []
        cellAreas = []
        min = 100
        for triangle in self.triangles:
            A = self.vertices[triangle[0]]
            B = self.vertices[triangle[1]]
            C = self.vertices[triangle[2]]

            barycenters.append([(A[0]+B[0]+C[0])/3,(A[1]+B[1]+C[1])/3])
            cellAreas.append(abs( (B[0]-A[0])*(C[1]-A[1]) - (B[1]-A[1])*(C[0]-A[0]))/2)
            if(cellAreas[-1] < min):
                min = cellAreas[-1]
            elif(cellAreas[-1] == 0):
                print(A,B,C)
                raise ValueError("Degenerated mesh : at least one of the triangle is of area 0.")
        self.cellAreas = np.array(cellAreas)
        self.barycenters = np.array(barycenters)
        return(min)

    def computeEdgeLength(self) -> None:
        """
        Compute the edgeLength array.

        Note:
            Required for:

                - LWRSolver
                - EikonalSolver
        """
        edgeLength = []
        for edge in self.edges:
            A = self.vertices[edge[0]]
            B = self.vertices[edge[1]]

            edgeLength.append(np.sqrt( (B[0]-A[0])*(B[0]-A[0]) + (B[1]-A[1])*(B[1]-A[1])))
        self.edgeLength = np.array(edgeLength)

    def getLimits(self):
        """
        Computes and returns the extremal coordinates of the vertices of the mesh.

        Returns:
            List[List[float]]: the extremal coordinates of the vertices of the mesh as a list of two two-elements lists i.e. [[x_min,x_max],[y_min,y_max]].
        """
        x_min, x_max = self.vertices[0][0],self.vertices[0][0]
        y_min, y_max = self.vertices[0][1],self.vertices[0][1]
        for P in self.vertices:
            if P[0] > x_max:
                x_max = P[0]
            if P[1] > y_max:
                y_max = P[1]
            if P[0] < x_min:
                x_min = P[0]
            if P[1] < y_min:
                y_min = P[1]
        return [[x_min,x_max],[y_min,y_max]]

    def saveToJson(self, filename):
        """
        Save the mesh object in a .json file.

        Args:
            filename (str): the name of the file where the mesh will be saved. The extension .json is not needed in the filename.

        Raises:
            ImportError: the json module must be installed.
        """
        if not json:
            raise ImportError("Module json not or wrongly installed. Needed for the json methods")

        MeshDico = {"type": "triangular mesh"}
        MeshDico["dx"] = self.dx
        MeshDico["vertices"] = self.vertices.tolist()
        MeshDico["edges"] = self.edges.tolist()
        MeshDico["triangles"] = self.triangles.tolist()

        MeshDico["outerVertices"] = [self.vertices[i].tolist() for i in self.boundaryPoints]

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

        MeshDico["zones"] = self.zones

        with open(filename+"_mesh.json", 'w', encoding='utf-8') as f:
            json.dump(MeshDico, f, ensure_ascii=False, indent=4)

    def loadFromJson(self, filename):
        """
        Load the mesh object with the data from a .json file.

        Args:
            filename (str): the name of the file to load. The extension .json is needed in the filename.

        Raises:
            ImportError: the json module must be installed.
        """
        if not json:
            raise ImportError("Module json not or wrongly installed. Needed for the json methods")
        with open(filename) as f:
            data = json.load(f)
            self.dx : float = data['dx']
            print("Minimal area for a triangle in the mesh : ", self.dx)

            self.vertices : ArrayLike = np.array(data['vertices'])
            self.edges : ArrayLike = np.array(data['edges']) #Liste des edges par paires d'indices de vertex
            self.triangles : ArrayLike = np.array(data['triangles'])

            self.exitVertices : ArrayLike = np.array(data['exitVertices'])
            self.exitEdges : ArrayLike = np.array(data['exitEdges'])#Liste des edges d'exit par indice d'edge
            self.wallEdges : ArrayLike = np.array(data['wallEdges'])#Liste des wall edges par indice d'edge
            self.boundaryPoints : ArrayLike = np.array(data['boundaryPoints'])
            self.boundaryEdgesIndex : ArrayLike = np.array(data['boundaryEdgesIndex'])#Liste des indices des edges faisant le bord du domaine
            self.vertexFlags : list = data['vertexFlags']

            self.trianglesWithEdges : ArrayLike = np.array(data['trianglesWithEdges'])#Liste des triangles par indice d'edge
            self.pairsOfTriangles : List[List[int]] = data['pairsOfTriangles'] #Liste des triangles entourant les edges ordonnée comme la EdgeList par les indices des triangles.
            self.trianglesPerVertex : List[list] = data['trianglesPerVertex'] #Nested List ordered by the vertex indices with the following content : [[triangle index, [otherVertex1, otherVertex2 ]], ...]

            self.outerNormalVectByTriangles : ArrayLike = np.array(data['outerNormalVectByTriangles'])#Liste de triplets correspondants aux 3 vecteurs normaux unitaires dans le même ordre que les edges correspondant dans la TriangleWithEdgeList
            self.cellAreas : ArrayLike = np.array(data['cellAreas'])
            self.barycenters : ArrayLike = np.array(data['barycenters'])
            self.edgeLength : ArrayLike = np.array(data['edgeLength'])

            self.zones :dict = data["zones"]

            print("Mesh successfully loaded.")



    def show(self, with_domain: NonConvexDomain=None, preference = "plotly") -> None:
        """
        Plotting method for the Mesh object.

        Args:
            domain (NonConvexDomain, optional): domain object from which exits and walls can be recovered for plotting.
            preference (str, optional): set to "plotly" or "matplotlib" to chose the preferred plotting package. If only one package is installed the preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed and no preference is set, plotly is used.
        """
        if( go and (not plt or preference == "plotly") ):
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
        elif(plt):
            fig, ax = plt.subplots()
            self.addPlot(fig,ax)
            for edge in self.wallEdges:
                path = [    (self.vertices[self.edges[edge][0]][0], self.vertices[self.edges[edge][0]][1]),
                            (self.vertices[self.edges[edge][1]][0], self.vertices[self.edges[edge][1]][1]) ]
                ax.add_patch(patches.Polygon(path, edgecolor='black',linewidth=1))

            for exit in self.exitEdges:
                path = [    (self.vertices[self.edges[exit][0]][0], self.vertices[self.edges[exit][0]][1]),
                            (self.vertices[self.edges[exit][1]][0], self.vertices[self.edges[exit][1]][1]) ]
                ax.add_patch(patches.Polygon(path, edgecolor='red',linewidth=2))
            limits = self.getLimits()
            ax.set_xlim(limits[0][0],limits[0][1])
            ax.set_ylim(limits[1][0],limits[1][1])
            plt.axis('equal')
            plt.show()
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")

    def addPlot(self, fig, ax = None, preference = "plotly"):
        """
        Non-blocking plotting method for the mesh object.

        Args:
            fig (matplotlib.Figure or plotly.graph_objects.Figure): the instance of Figure to which the plot must be added. Depends on the library used.
            ax (matplotlib.Axes or None): the instance of Axes to use if matplotlib is used.
            preference (str): set to "plotly" or "matplotlib" to chose the preferred plotting package. If only one package is installed the preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed for this method to work. If both are installed and no preference is set, plotly is used.
        """
        if( go and (not plt or preference == "plotly") ):
            for T in self.triangles:
                fig.add_trace(go.Scatter(x=[self.vertices[i][0] for i in T]+[self.vertices[T[0]][0]],
                                        y=[self.vertices[i][1] for i in T]+[self.vertices[T[0]][1]],
                                fill="toself",
                                fillcolor="White",
                                mode="lines",
                                line=dict(
                                    color="Black",
                                    width=0.5
                                 )))
        elif(plt):
            for T in self.triangles:
                triangle = patches.Polygon([(self.vertices[i][0],self.vertices[i][1]) for i in T], edgecolor='black', facecolor='white')
                ax.add_patch(triangle)

        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")


class CellValueMap(object):
    """
    An object to represent a cell valued map that corresponds to a function which is constant on the triangles of a mesh (mostly densities).

    Attributes
    -----------

    Attributes:
        Mesh (Mesh): the corresponding Mesh object on which the map is defined.
        values (List[float]): the values on each triangle as a list ordered in the same way as Mesh.triangles.

    Methods
    ---------
    """

    def __init__(self, Mesh: Mesh):
        self.Mesh : Mesh = Mesh
        self.values : List[float] = [0 for _ in self.Mesh.triangles]

    def generateRandom(self,variability :float = 0.4, mean : float= 0.23):
        """
        Generates random values for the map on each triangle of the mesh.
        The random value follows a uniform law on [mean -variability/2, mean + variability/2].

        Args:
            variability (float): the length of the interval of random values.
            mean (float): the center of the interval of random values.
        """
        self.values = [ mean + variability*(alea.random()-0.5) for _ in self.Mesh.triangles]

    def __len__(self):
        return len(self.values)

    def __getitem__(self, index):
        return self.values[index]

    def __setitem__(self, index, value):
        self.values[index] = value

    def __add__(self, other):
        if isinstance(other, CellValueMap):
            result = CellValueMap(self.Mesh)
            for index in range(len(self.Mesh.triangles)):
                result[index] = self[index] + other[index]
            return result
        elif isinstance(other, float) or isinstance(other, int):
            result = CellValueMap(self.Mesh)
            for index in range(len(self.Mesh.triangles)):
                result[index] = self[index] + other
            return result
        else:
            raise TypeError("A CellValueMap can be summed only with another CellValueMap or a float or an integer.")

    def __mul__(self, other):
        if isinstance(other, CellValueMap):
            result = CellValueMap(self.Mesh)
            for index in range(len(self.Mesh.triangles)):
                result[index] = self[index]*other[index]
            return result
        elif isinstance(other, float) or isinstance(other, int):
            result = CellValueMap(self.Mesh)
            for index in range(len(self.Mesh.triangles)):
                result[index] = self[index]*other
            return result
        else:
            raise TypeError("A CellValueMap can only be multiplied with another CellValueMap or a float or an integer.")

    def integrate(self) -> float:
        """
        Computes the integral of the map over the whole domain.

        Returns:
            float: the integral of the CellValueMap over the whole domain.
        """
        return sum([self.values[i]*self.Mesh.cellAreas[i] for i in range(len(self.Mesh.triangles))])

    def setConstantCircle(self, center:List[PointType], radius:float, value:float) -> None:
        """
        Sets to the given value all the cells for which the barycenter is inside the disk of given center and radius.

        Args:
            center (List[PointType]): the center of the disk
            radius (float): the radius of the disk
            value (float): the value to set on the corresponding cells.
        """
        for index in range(len(self.Mesh.triangles)):
            if((self.Mesh.barycenters[index][0] - center[0])**2 + (self.Mesh.barycenters[index][1] - center[1])**2 <= radius**2):
                self.values[index] = value

    def convolutionOverSquareBall(self, radius:float, conv_func) -> list:
        """
        Computes the convolution of the CellValueMap with the conv_func function on the support defined by the square of given radius (infinity-norm ball).

        Args:
            radius (float): the radius of the convolution support.
            conv_func (function: float -> float): the convolution function in the form  F(rho(y),|x_1-y_1|,|x_2,y_2|) where rho is the cellValueMap and x = (x_1,x_2) is the vertex where the convolution is computed. The computed quantity is then, for any x = (x_1,x_2), iint_{[-radius/2, radius/2]^2} F(rho(x+y), |y_1|,|y_2|) d y_1 d y_2.

        Returns:
            list: a list containing the computed value of the convolution for each vertex of the mesh ordered in the same way as the vertices of the Mesh object.
        """
        def recursiveIntegral(center, rad, visited, index) -> float:
            visited.append(index)
            distX = np.abs(self.Mesh.barycenters[index][0] - center[0])
            distY = np.abs(self.Mesh.barycenters[index][1] - center[1])
            if(distX > rad or distY > rad):
                return 0.0
            Sum = conv_func(self.values[index],distX,distY)*self.Mesh.cellAreas[index]
            for vertex in self.Mesh.triangles[index]:
                for neighborTriangle in self.Mesh.trianglesPerVertex[vertex]:
                    if neighborTriangle[0] not in visited:
                        Sum += recursiveIntegral(center,rad, visited, neighborTriangle[0])
            return Sum

        return [recursiveIntegral(self.Mesh.vertices[vertexIndex],radius,[],self.Mesh.trianglesPerVertex[vertexIndex][0][0]) for vertexIndex in range(len(self.Mesh.vertices))]

    def fitAveragedMap(self, other):
        """
        Fits a cellValueMap over another cellValueMap when the Mesh are different but the domains are the same by averaging over each triangle.

        Args:
            other (CellValueMap): the other map to fit onto the self map.
        """
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

    def show(self, preference = "plotly"):
        """
        Plotting method for the CellValueMap object.

        Args:
            preference (str, optional): set to "plotly" or "matplotlib" to chose the preferred plotting package. If only one package is installed the preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed and no preference is set, plotly is used.

        """
        if( go and (not plt or preference == "plotly") ):
            fig = go.Figure()
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
        elif(plt):
            fig, ax = plt.subplots()

            col = collections.PolyCollection([[ (self.Mesh.vertices[i][0],self.Mesh.vertices[i][1]) for i in T] for T in self.Mesh.triangles])
            col.set_cmap(cm.viridis)
            col.set_clim([0, 1])
            rgcol = col.set_array(self.values)

            ax.add_collection(col)
            fig.colorbar(rgcol, ax=ax, label="density")

            limits = self.Mesh.getLimits()
            ax.set_xlim(limits[0][0],limits[0][1])
            ax.set_ylim(limits[1][0],limits[1][1])
            plt.axis('equal')
            plt.show()
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")

    def getScatter(self, preference="plotly"):
        """
        Non-blocking plotting method for the cellValueMap object.

        Args:
            preference (str, optional): set to "plotly" or "matplotlib" to chose the preferred plotting package. If only one package is installed the preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed and no preference is set, plotly is used.
        """
        if( go and (not plt or preference == "plotly") ):
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
        elif(plt):
            print("Useless method for matplotlib.")
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")


class VertexValueMap(object):
    """
    An object to represent a vertex valued map that corresponds to a function which is affine on the triangles, defined by its values on the vertices (often a potential).

    Attributes
    -----------

    Attributes:
        Mesh (Mesh): the corresponding Mesh object on which the map is defined.
        values (List[float]): the values on each vertex as a list ordered in the same way as Mesh.vertices.

    Methods
    ---------
    """

    def __init__(self, Mesh):
        self.Mesh = Mesh
        self.values = np.array([0 for _ in self.Mesh.vertices])

    def generateRandom(self,variability = 0.5):
        """
        Generates random values for the maps distributed as 0.5 + variability*X where X is uniform on [-0.5,0.5].

        Args:
            variability (float, optional): the range of the uniform distribution.
        """
        self.values = np.array([ 0.5 + variability*(alea.random()-0.5) for _ in self.Mesh.vertices])

    def __len__(self):
        return len(self.values)

    def __getitem__(self, index): #renvoie le tuple (triangle, value)
        return self.values[index]

    def __setitem__(self,index, value):
        self.values[index] = value

    def __add__(self,other):
        newOne = VertexValueMap(self.Mesh)
        newOne.values = self.values + other.values

    def __str__(self):
        return(str(self.values))

    def setInfinity(self):
        """
        Sets all values of the map to *float('inf')*.
        """
        self.values = [float('inf') for _ in self.Mesh.vertices]


    def show3D(self, preference="plotly"):
        """
        Displays a 3D plot of the vertex valued map.

        Note:
            Only available with plotly at the moment.

        Args:
            preference (str, optional): set to "plotly" or "matplotlib" to chose the preferred plotting package. If only one package is installed the preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed and no preference is set, plotly is used.
        """
        if( go and (not plt or preference == "plotly") ):
            fig = go.Figure()
            fig.add_trace(go.Mesh3d(x=[P[0] for P in self.Mesh.vertices],
                                    y = [P[1] for P in self.Mesh.vertices],
                                    z = self.values,
                                    opacity=1,
                                    color='rgba(244,22,100,0.6)'))
            fig.show()
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")

    def add3Dplot(self, fig, color=[244,22,100,0.6], preference ="plotly"):
        """
        Non-blocking method that adds the 3D plot to a given `Figure` object.

        Note:
            Only available with plotly at the moment.

        Args:
            fig (matplotlib.pyplot.Figure or plotly.graph_objects.Figure): the instance of Figure to which the plot must be added. Depends on the library used.
            color (List[float]): color of the plot, the RGBA code.
            preference (str, optional): set to "plotly" or "matplotlib" to chose the preferred plotting package. If only one package is installed the preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed and no preference is set, plotly is used.

        """
        if( go and (not plt or preference == "plotly") ):
            fig.add_trace(go.Mesh3d(x=[P[0] for P in self.Mesh.vertices],
                                    y = [P[1] for P in self.Mesh.vertices],
                                    z = self.values,
                                    opacity=1,
                                    color='rgba('+str(color[0])+','+str(color[1])+','+str(color[2])+','+str(color[3])+')'))
        else:
            raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")


    def show(self, grid=False, colorscale_name = 'viridis', preference="plotly"):
        """
        Displays the vertex valued map as a colorscale scatter plot in 2D.

        Note:
            Only available with plotly at the moment.

        Args:
            grid (bool, optional): determines if the mesh is plotted or not.
            colorscale_name (str, optional): name of the colorscale to use for the plot.
            preference (str, optional): set to "plotly" or "matplotlib" to chose the preferred plotting package. If only one package is installed the preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed and no preference is set, plotly is used.

        """
        if( go and (not plt or preference == "plotly") ):
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


    def computeGradientFlow(self,normalize:bool = True, normalization = (lambda x,y : np.sqrt(x**2 + y**2))):
        """
        Computes and returns the gradients of the function defined by the vertex valued map. The gradient is computed on each triangle.

        Args:
            normalize (bool, optional): determines if the gradients should be renormalized.
            normalization (function: float, float -> float): the norm to use for the renormalization if `normalize=True`.

        Raises:
            ValueError: if at least one of the triangles is degenerated.

        Returns:
            ArrayLike: an array containing all the gradients ordered as the list Mesh.triangles.
        """
        LTrianglesGrad = []

        for triangle in self.Mesh.triangles: #On a des déterminants égaux à 0....
                det =  ((self.Mesh.vertices[triangle[1]][0] - self.Mesh.vertices[triangle[0]][0])*(self.Mesh.vertices[triangle[2]][1] - self.Mesh.vertices[triangle[0]][1])
                        - (self.Mesh.vertices[triangle[2]][0] - self.Mesh.vertices[triangle[0]][0])*(self.Mesh.vertices[triangle[1]][1] - self.Mesh.vertices[triangle[0]][1]))
                if(det == 0):
                    raise ValueError("At least one of the triangles of the mesh is degenerated.")
                Vecx = ( (self.values[triangle[0]] - self.values[triangle[2]])*(self.Mesh.vertices[triangle[1]][1] - self.Mesh.vertices[triangle[0]][1])
                        + (self.values[triangle[1]] - self.values[triangle[0]])*(self.Mesh.vertices[triangle[2]][1] - self.Mesh.vertices[triangle[0]][1]) )/det
                Vecy = ( (self.values[triangle[0]] - self.values[triangle[1]])*(self.Mesh.vertices[triangle[2]][0] - self.Mesh.vertices[triangle[0]][0])
                        + (self.values[triangle[2]] - self.values[triangle[0]])*(self.Mesh.vertices[triangle[1]][0] - self.Mesh.vertices[triangle[0]][0]) )/det
                LTrianglesGrad.append([-Vecx/normalization(Vecx, Vecy),-Vecy/normalization(Vecx, Vecy)] if normalize else [-Vecx,-Vecy])
        return np.array(LTrianglesGrad)

    def computeVertexGradientFlow(self,normalize:bool = True, normalization = (lambda x,y : np.sqrt(x**2 + y**2))):
        """
        Computes and returns the gradients of the function defined by the vertex valued map. The gradient is computed at each vertex as a weighted mean of the differences with the neighbouring vertices.

        Note:
            This method is less precise than `computeGradientFlow` in general. However due to its non-local construction, the gradient flow obtained with this method tends to be more regular.

        Args:
            normalize (bool, optional): determines if the gradients should be renormalized.
            normalization (function: float, float -> float): the norm to use for the renormalization if `normalize=True`.

        Returns:
            ArrayLike: an array containing all the gradients ordered as the list Mesh.triangles.
        """
        VertexTrianglesGrad = []

        for vertex in range(len(self.Mesh.vertices)):
            meanVect = np.array([0.0,0.0])
            treatedVertices = []
            for i in range(len(self.Mesh.trianglesPerVertex[vertex])):
                for point in self.Mesh.trianglesPerVertex[vertex][i][1]:
                    if point not in treatedVertices:
                        treatedVertices.append(point)
                        meanVect += ((self.values[point] - self.values[vertex])/np.sqrt((self.Mesh.vertices[point][0]-self.Mesh.vertices[vertex][0])**2 + (self.Mesh.vertices[point][1]-self.Mesh.vertices[vertex][1])**2))*(self.Mesh.vertices[point]-self.Mesh.vertices[vertex])
            VertexTrianglesGrad.append([meanVect[0]/normalization(meanVect[0], meanVect[1]),meanVect[1]/normalization(meanVect[0], meanVect[1])] if normalize else meanVect/len(treatedVertices))
        MeanTriangleGrads = []
        for triangle in self.Mesh.triangles:
            Vecx = sum([VertexTrianglesGrad[i][0] for i in triangle])
            Vecy = sum([VertexTrianglesGrad[i][1] for i in triangle])
            MeanTriangleGrads.append([-Vecx/normalization(Vecx, Vecy),-Vecy/normalization(Vecx, Vecy)] if normalize else [-Vecx,-Vecy])
        return np.array(MeanTriangleGrads)


    def showVectorField(self, normalize:bool = True, normalization = (lambda x,y : np.sqrt(x**2 + y**2)), preference="plotly"):
        """
        Displays the vector field corresponding to the gradient flow (obtained with `computeGradientFlow`) of the vertex valued map.

        Note:
            Only available with plotly at the moment.

        Args:
            normalize (bool, optional): determines if the gradients should be renormalized.
            normalization (function: float, float -> float): the norm to use for the renormalization if `normalize=True`.
            preference (str, optional): set to "plotly" or "matplotlib" to chose the preferred plotting package. If only one package is installed the preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed and no preference is set, plotly is used.
        """
        if( go and (not plt or preference == "plotly") ):
            L = self.computeGradientFlow(normalize, normalization)
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
