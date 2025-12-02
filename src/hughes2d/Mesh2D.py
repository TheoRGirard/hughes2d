"""
© Copyright 2025 Girard Théo

This file is part of the Hughes2d package.

The Hughes2d package is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.

The Hughes2d package is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the
Hughes2d package. If not, see <https://www.gnu.org/licenses/>.
"""

import random as alea
from pathlib import Path

import numpy as np
import triangle as tr

FigureType = None
AxesType = None

try:
    import matplotlib.pyplot as plt
    from matplotlib import cm, collections, patches
    #FigureType = FigureType | plt.Figure
    #AxesType = AxesType | plt.Axes
except ImportError:
    plt = None

try:
    import plotly.figure_factory as ff
    import plotly.graph_objects as go
    #FigureType = FigureType | go.Figure
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

#for Type checking
from collections.abc import Callable
from typing import Self

from numpy.typing import ArrayLike

PointType = list[float]

PRECISION = 1e-10
DEFAULT_MESH_DICT = {"domain" : 0, "exit" : [98], "wall" : [99]}
PLACEHOLDER_DOMAIN = [[0,0],[1,0],[0,1]]
DEFAULT_REQUIREMENTS = ["all"]
DEFAULT_COLOR_RGBA = [244,22,100,0.6]
TWO = 2


class NonConvexDomain:
    """Object corresponding to the geometry of a domain. The domain can have walls and exits on its boundary as well as walls inside the domain.

    Examples:
        Creating a square domain::

            MyDomain = NonConvexDomain([[0,0],[0,3],[3,3],[3,0]])

        Adding an exit::

            MyDomain.add_exit([[0,1],[0,2]])

        Creating a hole in the center of the domain::

            MyDomain.add_wall([[1,1],[1,2],[2,2],[2,1]], cycle=True)

    Args:
        outer_vertices_list (list[list[float]): The ordered list of the coordinates of
            the vertices defining the boundary of the domain. A vertex is represented
            as a pair of floats i.e. [float,float].

    Attributes
    -----------

    Attributes:
        outer_vertices (list[list[float]]): The list of the coordinates of the vertices
            located on the boundary of the domain. A vertex is represented as a pair of
            floats i.e. [float,float].
        wall_vertices (list[list[float]]): The list of the coordinates of the vertices
            located on walls inside the domain. A vertex is represented as a pair of
            floats i.e. [float,float].
        outer_boundary (list[list[int]]): The list of the edges of the outer boundary
            of the domain. An edge is represented as a pair of indices of the vertices
            of outer_vertices i.e. [int,int].
        wall_edges (list[list[int]]): The list of the edges of the outer boundary of
            the domain. An edge is represented as a pair of indices of the vertices of
            outer_vertices i.e. [int,int].
        wall_holes_point (list[list[float]]): The list of the coordinates of points
            located inside the holes of the domains. A point is represented as a pair
            of floats i.e. [float,float].
        exit_list (list[int]): The list of the edges corresponding to exits of the
            domain. An edge is represented as a pair of indices of the vertices of
            outer_vertices i.e. [int,int].
        zones (dict): dictionary containing all the zones defined on the mesh.

    Methods
    --------
    """

    def __init__(self, outer_vertices_list:list[PointType] = PLACEHOLDER_DOMAIN) -> None:
        self.outer_vertices = outer_vertices_list
        self.outer_boundary = []
        for i in range(len(outer_vertices_list)):
            self.outer_boundary.append([i,(i+1)%(len(outer_vertices_list))])
        self.wall_vertices = []
        self.wall_edges = []
        self.wall_holes_point = []
        self.exit_list = []
        self.zones = {}

    def import_from_dxf(self, filename: str) -> int:
        """Import from a dxf file.

        The dxf file can contain three layers:

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
        if not ezdxf:
            msg = "ezdxf not or wrongly installed."
            raise ImportError(msg)
        doc = ezdxf.readfile(filename)

        self.outer_vertices = []
        self.outer_boundary = []
        self.wall_vertices = []
        self.wall_edges = []
        self.wall_holes_point = []

        # iterate over all entities in modelspace
        msp = doc.modelspace()

        # entity query for all LINE entities in modelspace
        for e in msp.query("LINE"):
            if(e.dxf.layer not in {"0","domain"}):
                continue
            point1 = [np.round(e.dxf.start[0],5),np.round(e.dxf.start[1],5)]
            if point1 not in self.outer_vertices:
                id_point1 = len(self.outer_vertices)
                self.outer_vertices.append(point1)
            else:
                id_point1 = self.add_boundary_point(point1)

            point2 = [np.round(e.dxf.end[0],5),np.round(e.dxf.end[1],5)]
            if point2 not in self.outer_vertices:
                id_point2 = len(self.outer_vertices)
                self.outer_vertices.append(point2)
            else:
                id_point2 = self.add_boundary_point(point2)

            self.outer_boundary.append([id_point1, id_point2])

        if([len(self.outer_vertices)-1, 0] not in self.outer_boundary ):
            msg = "Corrupted domain : domain not closed."
            raise ValueError(msg)

        wall_vert = []
        for e in msp.query("LINE"):
            if(e.dxf.layer != "innerWalls"):
                continue
            cycling = 0
            point1 = [np.round(e.dxf.start[0],5),np.round(e.dxf.start[1],5)]
            if point1 not in wall_vert:
                wall_vert.append(point1)
            else:
                cycling += 1

            point2 = [np.round(e.dxf.end[0],5),np.round(e.dxf.end[1],5)]
            if point2 not in wall_vert:
                wall_vert.append(point2)
            else:
                cycling += 1

            if(cycling == 0):
                #New cycle !
                wall_vert.remove(point1)
                wall_vert.remove(point2)
                if(len(wall_vert) > 0):
                    self.add_wall(wall_vert, cycle =False)
                wall_vert = [point1,point2]
            elif(cycling == TWO):
                self.add_wall(wall_vert, cycle =True)
                wall_vert = []

        if(len(wall_vert) > 0):
            self.add_wall(wall_vert, cycle =False)

        for e in msp.query("LINE"):
            if(e.dxf.layer != "exits"):
                continue
            point1 = [np.round(e.dxf.start[0],5),np.round(e.dxf.start[1],5)]
            point2 = [np.round(e.dxf.end[0],5),np.round(e.dxf.end[1],5)]
            self.add_exit([point1,point2])

        zone_dict = {}
        for e in msp.query("LINE"):
            if(str(e.dxf.layer).find("zone_") != 0):
                continue
            label = str(e.dxf.layer)[5:]
            if(label not in zone_dict):
                zone_dict[label] = []

            point1 = [np.round(e.dxf.start[0],5),np.round(e.dxf.start[1],5)]
            if point1 not in zone_dict[label]:
                zone_dict[label].append(point1)

            point2 = [np.round(e.dxf.end[0],5),np.round(e.dxf.end[1],5)]
            if point2 not in zone_dict[label]:
                zone_dict[label].append(point2)

        for label, zone_vertices in zone_dict.items():
            self.add_zone(label, zone_vertices)



    def __contains__(self, point: PointType) -> bool:
        """Test if the point passed as a parameter is inside the convex hull of the domain.

        Args:
            point (list[float]): The coordinates of the point to test respresented
                as [float, float].

        Return:
            bool

        """
        for i in range(1,len(self.outer_vertices) - 1):
            if(belong_triangle(point, [  self.outer_vertices[0],
                                        self.outer_vertices[i],
                                        self.outer_vertices[i+1]])):
                return True
        return False

    def add_boundary_point(self, point : PointType) -> int:
        """Add a given point of the boundary to the list of points of the domain.

        This is typically used in order to guarantee that the given point will be
        included as a vertex of the mesh generated from this domain.
        If the point is already in outer_vertices, the method returns the index of
        the point without adding it.

        Args:
            point (list[float]): The point to add to the boundary respresented
                as [float, float].

        Raises:
            ValueError: The method raises a value error if the given point is not in
                the boundary of the domain.

        Returns:
            int: The method returns the index corresponding to the added point in the
                outer_vertices list.

        """
        if point not in self.outer_vertices:
            for index, edge in enumerate(self.outer_boundary):
                if(NonConvexDomain.belong_segment(point,
                                     [self.outer_vertices[edge[i]] for i in [0,1]])):
                    id_point = len(self.outer_vertices)
                    self.outer_vertices.append(point)
                    self.outer_boundary = [*(self.outer_boundary[:index]),
                                           [edge[0],id_point],
                                           [id_point, edge[1]],
                                           *(self.outer_boundary[(index+1):])]
                    return id_point
            msg = "The point "+ str(point) + " is not in the outer boundary."
            raise ValueError(msg)
        for i, point2 in enumerate(self.outer_vertices):
            if point2 == point:
                return i
        msg =   ("_point " + str(point) + " detected inside the boundary but doesn't"
                "correspond to any of the existing vertices.")
        raise ValueError(msg)

    def get_limits(self) -> list[PointType]:
        """Compute and return the extremal bounds of the domain.

        Returns:
            list[list[float]]: the extremal coordinates of the domain as a list of two
                two-elements lists i.e. [[x_min,x_max],[y_min,y_max]].

        """
        x_min, x_max = self.outer_vertices[0][0],self.outer_vertices[0][0]
        y_min, y_max = self.outer_vertices[0][1],self.outer_vertices[0][1]
        for point in self.outer_vertices:
            x_max = max(x_max,point[0])
            y_max = max(y_max,point[1])
            x_min = min(x_min,point[0])
            y_min = min(y_min,point[1])
        return [[x_min,x_max],[y_min,y_max]]


    def add_zone(self, zone_name:str, zone_vertices:list[PointType]) -> None:
        """Add a zone to the dict of zones of the domain.

        Args:
            zone_name (str): the name to be given to the new zone.
            zone_vertices (list[list[float]]): list of the vertex of the boundary of
                the zone.

        Raises:
            NameError: raises a name error if the zone name is already taken.

        """
        if(zone_name in self.zones):
            msg = "The zone name " + zone_name + " is already used."
            raise NameError(msg)
        self.zones[zone_name] = zone_vertices


    def add_wall_point(self, point: PointType) -> int:
        """Add a given point of the boundary to the list of points of the inner walls of the domain.

        This is typically used in order to guarantee that the given point will be
        included as a vertex of the mesh generated from this domain.
        If the point is already in wall_vertices, the method returns the index of the
        point without adding it.

        Args:
            point (list[float]): the point to add on an inner wall of the domain.

        Raises:
            ValueError: the method raises a value error if the given point is not in
                the inner walls of the domain.

        Returns:
            int: the method returns the index corresponding to the added point in the
                wall_vertices list.

        """
        if point not in self.wall_vertices:
            for index, edge in enumerate(self.wall_edges):
                if(NonConvexDomain.belong_segment(point,
                                     [self.wall_vertices[edge[i]] for i in [0,1]])):
                    id_point = len(self.wall_vertices)
                    self.wall_vertices.append(point)
                    self.wall_edges = [*self.wall_edges[:index],
                                       [edge[0],id_point],
                                       [id_point, edge[1]],
                                       *self.wall_edges[(index+1):]]
                    return id_point
            msg = "The point " + str(point) + " given is not in an inner wall."
            raise ValueError(msg)
        for i, point2 in enumerate(self.wall_vertices):
            if point2 == point:
                return i
        msg =   ("_point " + str(point) + " detected inside an inner wall but doesn't"
                "correspond to any of the existing vertices.")
        raise ValueError(msg)

    def add_wall(self, coord_wall: list[PointType],*, cycle: bool=False ) -> None:
        """Add a wall in the domain.

        If cycle is set to True, the wall is considered as an area to exclude from
        the domain. If not, only the edges defined by the coord_wall are excluded
        from the domain.

        Note:
            The cycle=True works properly only with convex holes in the domain. In fact,
                the barycenters of all the points of the walls to add must be inside
                the hole for the exclusion to work properly.

        Args:
            coord_wall (list[PointType]) : list of the vertices defining the walls.
            cycle (bool, keyword-only) : A boolean switching between a hole to exclude
                from the domain and only walls of zero thickness to exclude from
                the domain.

        Raises:
            ValueError: raises a Value Error if the points of wall_coord are not inside
                the convex hull of the domain.

        """
        for point in coord_wall:
            if(point not in self):
                msg = "The wall " + str(coord_wall) + " is not inside the domain."
                raise ValueError(msg)

        if(cycle): #CircularWall
            for i in range(len(self.wall_vertices),
                           len(self.wall_vertices)+len(coord_wall)-1):
                self.wall_edges.append([i,i+1])
            self.wall_edges.append([len(self.wall_vertices)+len(coord_wall)-1,
                                    len(self.wall_vertices)])

            for point in coord_wall:
                self.wall_vertices.append(point)

            center_point = [0,0]
            for point in coord_wall:
                center_point[0] += point[0]
                center_point[1] += point[1]

            self.wall_holes_point.append([center_point[0]/len(coord_wall),center_point[1]/len(coord_wall)])
        else:
            for i in range(len(self.wall_vertices),
                           len(self.wall_vertices)+len(coord_wall)-1):
                self.wall_edges.append([i,i+1])

            for point in coord_wall:
                self.wall_vertices.append(point)

    def add_exit(self, exit_edge: list[PointType]) -> None:
        """Add an exit passed in parameter.

        If the two exit points are on different edges, the shortest path (in number of
        vertices crossed) following the wall is added as an exit.

        Args:
            exit_edge (list[PointType]) : an exit defined by two extremal points.

        Raises:
            ValueError: raises an error if:

                - at least one of the extremal points is not in a wall edge or a
                    boundary edge;
                - there exists no path staying either in the walls or in the boundary
                    that links the two extremal points of the exit.

        """
        vertex_index = [0,0]
        if(len(exit_edge) > 2):
            msg = "Invalid value: "+str(exit_edge)+", an exit is a list of two points."
            raise ValueError(msg)
        if(self.has_wall_edge(exit_edge)):
            for i in [0,1]:
                vertex_index[i] = self.add_wall_point(exit_edge[i])

            exit_path = NonConvexDomain.shortest_path_bfs(vertex_index[0],
                                                        vertex_index[1],
                                                        self.wall_edges)
            if len(exit_path) == 0:
                msg = "Impossible exit for this domain inner walls"
                raise ValueError(msg)

            exit_edges = [ [self.wall_vertices[exit_path[i]],
                            self.wall_vertices[exit_path[i+1]]]
                          for i in range(len(exit_path) - 1)]
            self.exit_list += exit_edges

        elif(self.has_outer_edge(exit_edge)):
            for i in [0,1]:
                vertex_index[i] = self.add_boundary_point(exit_edge[i])
            exit_path = NonConvexDomain.shortest_path_bfs(vertex_index[0],
                                                        vertex_index[1],
                                                        self.outer_boundary)
            if len(exit_path) == 0:
                msg = "Impossible exit for this domain boundary"
                raise ValueError(msg)

            exit_edges = [ [self.outer_vertices[exit_path[i]],
                            self.outer_vertices[exit_path[i+1]]]
                          for i in range(len(exit_path) - 1)]
            self.exit_list += exit_edges
        else:
            msg = ("The exit "+str(exit_edge)+" is not inside an inner wall nor in"
                    "an outer oboundary.")
            raise ValueError(msg)

    def add_exits(self, list_of_exits: list[list[PointType]]) -> None:
        """Call the method add_exit multiple times.

        Args:
            list_of_exits (list[list[PointType]]): a list of the exits to add.

        """
        for exit_edge in list_of_exits:
            self.add_exit(exit_edge)

    def find_boundary_point(self, point: PointType) -> int :
        """Search an edge of the boundary containing the given point.

        Args:
            point (PointType): the point the should be in the researched edge.

        Returns:
            int: the index of the first edge found in the boundary that contains P;
                returns -1 if no such edge is found.

        """
        for i,edge in enumerate(self.outer_boundary):
            if NonConvexDomain.belong_segment(point, [self.outer_vertices[edge[0]],
                                                     self.outer_vertices[edge[1]]]):
                return i
        return -1

    def find_wall_point(self, point: PointType) -> int :
        """Search an edge of the inner walls containing a given point.

        Args:
            point (list[float]): the point the should be in the researched edge.

        Returns:
            int: the index of the first edge found in the inner walls that contains
                the point; returns -1 if no such edge is found.

        """
        for i, edge in enumerate(self.wall_edges):
            if NonConvexDomain.belong_segment(point, [self.wall_vertices[edge[0]],
                                                     self.wall_vertices[edge[1]]]):
                return i
        return -1

    def find_exit_point(self, point: PointType) -> int:
        """Search an edge of the exits containing a given point.

        Args:
            point (list[float]): the point the should be in the researched edge.

        Returns:
            int: the index of the first edge found in the exits that contains the point;
                returns -1 if no such edge is found.

        """
        for i, exit_edge in enumerate(self.exit_list):
            if NonConvexDomain.belong_segment(point, exit_edge):
                return i
        return -1

    def has_boundary_point(self, point: PointType) -> bool :
        """Check that the point passed as a parameter belongs to the outer boundary.

        Args:
            point (list[float]): the point to test.

        Returns:
            bool: A boolean that is True if the point is in the outer boundary.

        """
        return (self.find_boundary_point(point) != -1)

    def has_wall_point(self, point: PointType) -> bool :
        """Check that the point passed as a parameter belongs to the inner walls.

        Args:
            point (list[float]): the point to test.

        Returns:
            bool: A boolean that is True if the point is in the inner walls.

        """
        return (self.find_wall_point(point) != -1)

    def has_exit_point(self, point: PointType) -> bool:
        """Check that the point passed as a parameter belongs to the exits.

        Args:
            point (list[float]): the point to test.

        Returns:
            bool: A boolean that is True if the point is in the exits.

        """
        return (self.find_exit_point(point) != -1)

    def has_exit_edge(self, edge : list[PointType]) -> bool:
        """Test if the segment defined by the pair of points passed as a parameter is inside an exit.

        Args:
            edge (list[PointType]) :  a pair of points defining a segment to test.

        Returns:
            bool

        """
        exit_indices = [self.find_exit_point(edge[0]), self.find_exit_point(edge[1])]
        if exit_indices[0] == -1 or exit_indices[1] == -1:
            return False
        if exit_indices[0] == exit_indices[1]:
            return True
        return(  (  NonConvexDomain.belong_segment(   self.exit_list[exit_indices[0]][0],
                                            self.exit_list[exit_indices[1]])
                or NonConvexDomain.belong_segment(   self.exit_list[exit_indices[0]][1],
                                             self.exit_list[exit_indices[1]]) )
            and np.abs( (self.exit_list[exit_indices[0]][1][0]
                         - self.exit_list[exit_indices[0]][0][0])
                        *(self.exit_list[exit_indices[1]][1][1]
                          - self.exit_list[exit_indices[1]][0][1])
                        - (self.exit_list[exit_indices[0]][1][1]
                           - self.exit_list[exit_indices[0]][0][1])
                        *(self.exit_list[exit_indices[1]][1][0]
                          - self.exit_list[exit_indices[1]][0][0]) )
                  < PRECISION )

    def has_wall_edge(self, edge: list[PointType]) -> bool:
        """Test if the segment defined by the pair of points passed as a parameter is inside an inner wall.

        Args:
            edge (list[PointType]) :  a pair of points defining a segment to test.

        Returns:
            bool

        """
        wall_indices = [self.find_wall_point(edge[0]), self.find_wall_point(edge[1])]
        if wall_indices[0] == -1 or wall_indices[1] == -1:
            return False
        wall1_coord = [self.wall_vertices[self.wall_edges[wall_indices[0]][0]],
                       self.wall_vertices[self.wall_edges[wall_indices[0]][1]]]
        wall2_coord = [self.wall_vertices[self.wall_edges[wall_indices[1]][0]],
                       self.wall_vertices[self.wall_edges[wall_indices[1]][1]]]
        if (NonConvexDomain.belong_segment(edge[0],wall2_coord)
            or NonConvexDomain.belong_segment(edge[1],wall1_coord)):
            return True
        if (    (NonConvexDomain.belong_segment(wall1_coord[0], wall2_coord)
                 or NonConvexDomain.belong_segment(wall1_coord[1], wall2_coord))
            and np.abs( (wall1_coord[1][0]- wall1_coord[0][0])
                       *(wall2_coord[1][1]- wall2_coord[0][1])
                       - (wall1_coord[1][1]- wall1_coord[0][1])
                       *(wall2_coord[1][0]- wall2_coord[0][0]) )
               < PRECISION):
                return True

        vertices_index = [0,0]
        for index, point in enumerate(self.wall_vertices):
            for i in range(2):
                if (((edge[i][0] - point[0])**2 + (edge[i][1] - point[1])**2)
                        < PRECISION ):
                    vertices_index[i] = index
        for wall in self.wall_edges:
            if( (wall[0] == vertices_index[0] and wall[1] == vertices_index[1])
                or (wall[0] == vertices_index[1] and wall[1] == vertices_index[0])):
                return True

        return False

    def has_outer_edge(self, edge : list[PointType]) -> bool:
        """Test if the segment defined by the pair of points passed as a parameter is inside the outer boundary.

        Args:
            edge (list[PointType]) :  a pair of points defining a segment to test.

        Returns:
            bool

        """
        outer_indices = [self.find_boundary_point(edge[0]), self.find_boundary_point(edge[1])]
        if outer_indices[0] == -1 or outer_indices[1] == -1:
            return False
        outer1_coord = [self.outer_vertices[self.outer_boundary[outer_indices[0]][0]],
                        self.outer_vertices[self.outer_boundary[outer_indices[0]][1]]]
        outer2_coord = [self.outer_vertices[self.outer_boundary[outer_indices[1]][0]],
                        self.outer_vertices[self.outer_boundary[outer_indices[1]][1]]]
        if( NonConvexDomain.belong_segment(edge[0],outer2_coord)
           or NonConvexDomain.belong_segment(edge[1],outer1_coord)):
            return True
        if( ( NonConvexDomain.belong_segment(outer1_coord[0], outer2_coord)
             or NonConvexDomain.belong_segment(outer1_coord[1], outer2_coord))
            and np.abs( (outer1_coord[1][0]- outer1_coord[0][0])
                       *(outer2_coord[1][1]- outer2_coord[0][1])
                       - (outer1_coord[1][1]- outer1_coord[0][1])
                       *(outer2_coord[1][0]- outer2_coord[0][0]) )
                   < PRECISION):
            return True
        for i in self.outer_boundary[outer_indices[0]]:
            for j in self.outer_boundary[outer_indices[1]]:
                if (self.outer_vertices[i][0] == edge[0][0]
                     and self.outer_vertices[i][1] == edge[0][1]
                     and self.outer_vertices[j][0] == edge[1][0]
                     and self.outer_vertices[j][1] == edge[1][1]):
                     for edge_test in self.outer_boundary:
                            if ((edge_test[0] == i and edge_test[1] == j)
                                or (edge_test[1] == i and edge_test[0] == j)):
                                return True
        return False

    def show(self, preference:str ="plotly") -> None:
        """Plot method for the domain object.

        The method opens either a new window or a default navigator tab.

        Args:
            preference (str, optional): set to "plotly" or "matplotlib" to chose the
                preferred plotting package. If only one package is installed the
                preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed for this method
                to work. If both are installed and no preference is set, plotly is used.

        """
        if( go and (not plt or preference == "plotly") ):
            fig = go.Figure()
            self.add_plot(fig,ax=None,preference=preference)
            fig.update_layout(yaxis={
                "scaleanchor" : "x",
                "scaleratio" : 1})
            fig.show()
        elif(plt):
            fig, ax = plt.subplots()
            self.add_plot(fig,ax,preference)
            limits = self.get_limits()
            ax.set_xlim(limits[0][0],limits[0][1])
            ax.set_ylim(limits[1][0],limits[1][1])
            plt.axis("equal")
            plt.show()
            fig.clf()
            plt.close(fig)
        else:
            msg = ("No plotting module found. Try installing plotly or matplotlib"
                   " if you want to use show methods")
            raise ImportError(msg)

    def add_plot(self, fig: FigureType, ax: AxesType = None, preference: str = "plotly") -> None:
        """Non-blocking plot method for the domain object, does not show the graph.

        Args:
            fig (matplotlib.Figure or plotly.graph_objects.Figure): the instance of
                Figure to which the plot must be added. Depends on the library used.
            ax (matplotlib.Axes or None): the instance of Axes to use if matplotlib is
                used.
            preference (str): set to "plotly" or "matplotlib" to chose the preferred
                plotting package. If only one package is installed the preference is
                ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed for this method to
                work. If both are installed and no preference is set, plotly is used.

        """
        if( go and (not plt or preference == "plotly") ): #plotly version of the plot
            start_point = self.outer_vertices[self.outer_boundary[0][0]]
            ordered_outer_vertices = [start_point]
            end_point = self.outer_vertices[self.outer_boundary[0][1]]
            while(end_point != start_point):
                ordered_outer_vertices.append(end_point)
                for edge in self.outer_boundary:
                    if(self.outer_vertices[edge[0]] == end_point):
                        end_point = self.outer_vertices[edge[1]]
                        break

            fig.add_trace(go.Scatter(x=[ *[P[0] for P in ordered_outer_vertices],
                                        ordered_outer_vertices[0][0] ],
                                    y=[ *[P[1] for P in ordered_outer_vertices],
                                       ordered_outer_vertices[0][1]],
                                    fill="toself", fillcolor="White", mode="lines"))
            for edge in self.wall_edges:
                fig.add_shape(type="line",
                        x0=self.wall_vertices[edge[0]][0],
                        y0=self.wall_vertices[edge[0]][1],
                        x1=self.wall_vertices[edge[1]][0],
                        y1=self.wall_vertices[edge[1]][1],
                        line={
                            "color" : "LightSeaGreen",
                            "width" : 2,
                            })
            for exit_edge in self.exit_list:
                fig.add_shape(type="line",
                        x0=exit_edge[0][0],
                        y0=exit_edge[0][1],
                        x1=exit_edge[1][0],
                        y1=exit_edge[1][1],
                        line={
                            "color" : "Red",
                            "width" : 2,
                        })
        elif(plt): #matplotlib version of the plot
            start_point = self.outer_vertices[self.outer_boundary[0][0]]
            ordered_outer_vertices = [start_point]
            end_point = self.outer_vertices[self.outer_boundary[0][1]]
            while(end_point != start_point):
                ordered_outer_vertices.append(end_point)
                for edge in self.outer_boundary:
                    if(self.outer_vertices[edge[0]] == end_point):
                        end_point = self.outer_vertices[edge[1]]
                        break

            domain_polygon=patches.Polygon([(P[0],P[1])
                                            for P in ordered_outer_vertices],
                                             edgecolor="black", facecolor="white")
            ax.add_patch(domain_polygon)

            for edge in self.wall_edges:
                path = [(self.wall_vertices[edge[0]][0], self.wall_vertices[edge[0]][1]),
                        (self.wall_vertices[edge[1]][0], self.wall_vertices[edge[1]][1])]
                ax.add_patch(patches.Polygon(path, edgecolor="black",linewidth=1))

            for exit_edge in self.exit_list:
                path = [    (exit_edge[0][0], exit_edge[0][1]),
                            (exit_edge[1][0], exit_edge[1][1]) ]
                ax.add_patch(patches.Polygon(path, edgecolor="red",linewidth=2))

        else:
            msg = ("No plotting module found. Try installing plotly or matplotlib"
                   " if you want to use show methods")
            raise ImportError(msg)

    @staticmethod
    def shortest_path_bfs(i_start: int , i_end: int, list_edges: list[list[int]]) -> list[int]:
        """Compute the shortest path between i_start and i_end (in number of vertices) on the network defined by the integers as vertice and the edges of list_edges.

        The method uses a Breadth First Search algorithm with memory of the pathes
        explored.

        Args:
            i_start (int) : the starting vertex of the searched path.
            i_end (int) : the ending vertex of the searched path.
            list_edges (list[list[int]]) : a list of the edges of the network symbolized
                by a list of two integers.

        Returns:
            list[int]: The shortest path found as a list of integers i.e. the list of
                the successive vertices defining the path. If no path is found the
                method returns an empty list.

        """
        visited = []
        pathes = [[i_start]]
        while(len(pathes) > 0):
            new_pathes = []
            for path in pathes:
                if path[-1] == i_end:
                    return path
                visited.append(path[-1])

                for edge in list_edges:
                    if path[-1] == edge[0] and edge[1] not in visited:
                        new_pathes.append([*path,edge[1]])
                    if path[-1] == edge[1] and edge[0] not in visited:
                        new_pathes.append([*path,edge[0]])
            pathes = new_pathes

        return []

    @staticmethod
    def belong_segment(point: PointType, segment: list[PointType]) -> bool:
        """Check that a point belongs to a segment within an error margin of PRECISION.

        Args:
            point (PointType) : the point to test.
            segment (list[PointType]) : the segment to test.

        Returns:
            bool

        """
        point_a = segment[0]
        point_b = segment[1]
        if(np.abs((point_a[0]- point[0])*(point_a[1]-point_b[1])
                  -(point_a[1]- point[1])*(point_a[0]-point_b[0]))
           < PRECISION):
            scal =( (point_a[0] - point[0])*(point_a[0] - point_b[0])
                    +(point_a[1] - point[1])*(point_a[1] - point_b[1]) )
            if( scal <= (point_a[0] - point_b[0])**2 + (point_a[1] - point_b[1])**2
               and scal >= 0):
                return True
        return False

def belong_triangle(point: PointType,triangle: list[PointType]) -> bool:
    """Check that a point belongs to a triangle.

    Args:
        point (PointType) : the point to test.
        triangle (list[PointType]) : the triangle to test.

    Returns:
        bool

    """
    det = ( (triangle[1][0]-triangle[0][0])*(triangle[2][1]-triangle[0][1])
            - (triangle[1][1]-triangle[0][1])*(triangle[2][0]-triangle[0][0]) )
    if det == 0:
        return False

    projection_x = ((point[0]-triangle[0][0])*(triangle[2][1]-triangle[0][1])
                     - (point[1]-triangle[0][1])*(triangle[2][0]-triangle[0][0]))/det
    projection_y = ((triangle[1][0]-triangle[0][0])*(point[1]-triangle[0][1])
                     - (triangle[1][1]-triangle[0][1])*(point[0]-triangle[0][0]))/det
    return( projection_x >= 0 and projection_y >= 0 and projection_x+projection_y <= 1)



class Mesh:
    """The Mesh Object is an object representing a triangular mesh of a given grain for a given domain object.

    It contains all the useful lists of edges, triangles and vertices for optimized
    computations of the numerical scheme involved in the hughes2d package.

    Examples:
        The mesh can be generated by using a NonConvexDomain as the basis of
        computations::

            MyDomain = NonConvexDomain([[0,0],[0,1],[1,1],[1,0]])
            MyMesh = Mesh()
            MyMesh.generate_mesh_from_domain(MyDomain,0.1)

        It is also possible to import a mesh from a .msh file::

            MyMesh = Mesh()
            MyMesh.import_mesh_from_msh("filename.msh")

        See also *import_mesh_from_msh_free_fem* and *import_from_lists* for different
        import methods. The mesh can then be exported (see *export_mesh_msh* and
        *export_mesh_msh_free_fem*) or saved as a .json file::

            MyMesh.save_to_json("filename")

        And can be loaded from this .json file::

            MyMesh2 = Mesh()
            MyMesh2.load_from_json("filename.json")

    Attributes
    -----------

    Attributes:
        min_cell_area (float): minimal area of a triangle in the mesh.
        max_edge_length (float): maximal length of an edge in the mesh.
        vertices (ArrayLike): array of all the vertices triangle_edge_coordinates.
        edges (ArrayLike): array of all the edges as [vertexIndex, vertexIndex].
        triangles (ArrayLike): array of all the triangles
            as [vertexIndex, vertexIndex, vertexIndex].
        exit_vertices (ArrayLike): array of the vertices index that belong to an exit.
        exit_edges (ArrayLike): array of the edges index where the edge belongs to
            an exit.
        wall_edges (ArrayLike): array of the edges index where the edge belongs to
            a wall.
        boundary_points (ArrayLike): array of the indices of the vertices that belong
            to the outer boundary of the mesh.
        boundary_edges_index (ArrayLike): array of the edges index where the edge
            belongs to the outer boundary.
        triangles_with_edges (ArrayLike): array of triangles ordered as self.triangles,
            elements as [edge_index, edge_index, edge_index].
        pairs_of_triangles (list[list[int]]): nested list of length 1 or 2, ordered as
            self.edges, elements as [triangle_index, triangle_index] or [triangle_index]
        triangles_per_vertex (list[list]): nested list, ordered as self.vertices, an
            element is a list of (number of triangles containing the vertex) elements
            as [[triangle index, [otherVertex1, otherVertex2 ]], ...].
        outer_normal_vect_by_triangles (ArrayLike): array order as self.tringles
            containing the unit normal vectors corresponding to the three edges of
            each triangle. The normal vectors are directed towards the exterior of
            the triangle.
        cell_areas (ArrayLike): array ordered as self.triangles containing the area of
            each triangle of the mesh.
        edge_length (ArrayLike):array ordered as self.edges containing the length of
            each edge.
        zones (dict): dictionary containing all the zones defined on the mesh.


    Methods
    --------
    """

    def __init__(self) -> None:
        self.min_cell_area : float = None
        self.max_edge_length : float = None

        self.vertices : ArrayLike
        self.edges : ArrayLike
        self.triangles : ArrayLike

        self.exit_vertices : ArrayLike
        self.exit_edges : ArrayLike
        self.wall_edges : ArrayLike
        self.boundary_points : ArrayLike
        self.boundary_edges_index : ArrayLike
        self.vertex_flags : list = []

        self.triangles_with_edges : ArrayLike
        self.pairs_of_triangles : list[list[int]]
        self.triangles_per_vertex : list[list]

        self.outer_normal_vect_by_triangles : ArrayLike
        self.cell_areas : ArrayLike
        self.barycenters : ArrayLike
        self.edge_length : ArrayLike

        self.zones : dict = {}

    def import_mesh_from_msh(self, filename :str, requirements : list[str] = DEFAULT_REQUIREMENTS, *, verbose : bool = False) -> None:
        """Import the data from a .msh file into the Mesh object in the GMSH format.

        Args:
            filename (string) : the path to the file to import.
            verbose (bool, optional): displays information in the console during the
                import.
            requirements (list[str], optional): a list containing the computations that
                will be done using the mesh. The possible values are:

                - EikonalSolver : the mesh will be used in order to solve an eikonal
                    equation.
                - LWRSolver : the mesh will be used in order to solve a scalar
                    conservation law.
                - all : all the possible computations will be done
                - integrate : the mesh will be used in order to compute some
                    integrals over the domain.
                - FreeFEM : the mesh will be used in order to import/export mesh
                    file for FreeFEM.

        Raises:
            ImportError: requires the python library meshio.

        """
        if not meshio:
            msg = "meshio must be installed in order to use .msh related methods."
            raise ImportError(msg)
        file_mesh = meshio.read(filename)
        self.vertices = np.array([ [P[0],P[1]] for P in file_mesh.points ])
        self.triangles = file_mesh.cells_dict["triangle"]

        self.compute_edge_list()
        self.exit_edges = []
        self.wall_edges = []

        special_edges_index = []
        if("line" in file_mesh.cells_dict):
            for index, edge in enumerate(self.edges):
                for outer_edge in file_mesh.cells_dict["line"]:
                    if( (edge[0] == outer_edge[0] and edge[1] == outer_edge[1])
                       or (edge[0] == outer_edge[1] and edge[1] == outer_edge[0]) ):
                        special_edges_index.append(index)
                        break
            if(len(special_edges_index) != len(file_mesh.cells_dict["line"])):
                msg = ("Error importing " +str(filename)+ " : some line cells are not"
                       " edges from the mesh.")
                raise ImportError(msg)

        exit_edges = []
        wall_edges = []
        if("hughes2d:special" in file_mesh.cell_data):
            for mesh_index, edge_index in enumerate(special_edges_index):
                if file_mesh.cell_data["hughes2d"][0][mesh_index] == 1:
                    exit_edges.append(edge_index)
                elif file_mesh.cell_data["hughes2d"][0][mesh_index] == TWO:
                    wall_edges.append(edge_index)
                else:
                    print("Warning : some special edges are neither wall nor exit.")
            self.exit_edges = np.array(exit_edges)
            self.wall_edges = np.array(wall_edges)

        else :
            print("WARNING : Mesh without walls imported from msh. Setting all special"
                  " edges as walls.")
            wall_edges = [edge_index for edge_index in special_edges_index
                                        if edge_index not in exit_edges]

            self.wall_edges = np.array(wall_edges)

        if(len(self.exit_edges) == 0):
            print("WARNING : Mesh without exits imported from msh.")
        if(len(special_edges_index) == 0 and len(self.exit_edges) == 0
                                        and len(self.wall_edges) == 0):
            msg = "Error importing "+str(filename)+" : Mesh without walls neither exits"
            raise ImportError(msg)
        if(len(special_edges_index) != len(self.exit_edges) + len(self.wall_edges)):
            print("WARNING : Losing some special edges during the import of .msh file.")

        print(f"Mesh imported. Contains {len(self.triangles)} triangles.")
        self.computations(verbose=verbose, requirements=requirements)

    def import_mesh_from_msh_free_fem(self, filename : str, flag_dict : dict = DEFAULT_MESH_DICT, requirements : list[str] = DEFAULT_REQUIREMENTS, *, verbose : bool = False) -> None:
        """Import the data from a mesh file constructed in FreeFEM into the Mesh object.

        The specific structure (inner walls, exits...) can be specified by specifying
        different flags in FreeFEM.

        Args:
            filename (str) : the path to the file to import.
            flag_dict (dict) : a dictionary describing the specific translation of the
                FreeFEM flag number. Must contain the keys domain, exit and wall.
            verbose (bool, optional): displays information in the console during the
                import.
            requirements (list[str], optional): a list containing the computations that
                will be done using the mesh. The possible values are:

                - EikonalSolver : the mesh will be used in order to solve an eikonal
                    equation
                - LWRSolver : the mesh will be used in order to solve a scalar
                    conservation law
                - all : all the possible computations will be done
                - integrate : the mesh will be used in order to compute some
                    integrals over the domain
                - FreeFEM : the mesh will be used in order to import/export mesh
                    file for FreeFEM

        """
        with Path(filename).open("r") as file:
            line = file.readline().split()
            nb_vertices, nb_triangles, nb_spe_edges = int(line[0]), int(line[1]), int(line[2])

            vertices = []
            exit_vertices = []
            wall_vertices = []
            for i in range(nb_vertices):
                line = file.readline().split()
                vertices.append([float(line[0]), float(line[1])])
                if(int(line[2]) in flag_dict["exit"]):
                    exit_vertices.append(i)
                elif(int(line[2]) in flag_dict["wall"]):
                    wall_vertices.append(i)
            self.vertices = np.array(vertices)
            self.exit_vertices = np.array(exit_vertices)
            self.boundary_points = np.array(exit_vertices+wall_vertices)

            triangles = []
            for _ in range(nb_triangles):
                line = file.readline().split()
                triangles.append([int(line[0])-1, int(line[1])-1, int(line[2])-1])
            self.triangles = np.array(triangles)

            self.compute_edge_list()

            exit_edges = []
            wall_edges = []

            for _ in range(nb_spe_edges):
                line = file.readline().split()
                for index, edge in enumerate(self.edges):
                    if( (edge[0] == int(line[0])-1 and edge[1] == int(line[1])-1 )
                       or (edge[0] == int(line[1])-1  and edge[1] == int(line[0])-1 ) ):
                        if(int(line[2]) in flag_dict["exit"]):
                            exit_edges.append(index)
                        elif(int(line[2]) in flag_dict["wall"]):
                            wall_edges.append(index)
                        else :
                            print("Warning : some special edges are neither wall nor exit")
                        break

            self.exit_edges = np.array(exit_edges)
            self.wall_edges = np.array(wall_edges)


            if(len(self.exit_edges) == 0):
                print("WARNING : Mesh without exits imported from msh.")
            if(len(self.exit_edges) == 0 and len(self.wall_edges) == 0):
                msg = ("Error importing "+str(filename)+" : Mesh without"
                       " walls neither exits")
                raise ImportError(msg)
            if(nb_spe_edges != len(self.exit_edges) + len(self.wall_edges)):
                print("WARNING : Losing some special edges during the import of .msh file.")

            if verbose :
                print(f"Mesh imported. Contains {len(self.triangles)} triangles.")

            self.computations(verbose=verbose, requirements=requirements)

    def import_from_lists(self, vertices:list[PointType], triangles:list[list[int]], domain:NonConvexDomain, requirements : list[str] = DEFAULT_REQUIREMENTS, *, verbose :bool = False) -> None:
        """Import the lists passed as parameters into the Mesh object.

        Args:
            vertices (list[PointType]) : list of all the vertices coordinates.
            triangles (list[list[int]]) : list of the triangles symbolized by a triplet
                of the indices of the corresponding verttices in the vertices list.
            domain (NonConvexDomain) : a NonConvexDomain instance corresponding to the
                mesh imported.
            verbose (bool, optional): displaying progression in console.
            requirements (list[str], optional): a list containing the computations that
                will be done using the mesh. The possible values are:

                - EikonalSolver : the mesh will be used in order to solve an eikonal
                    equation
                - LWRSolver : the mesh will be used in order to solve a scalar
                    conservation law
                - all : all the possible computations will be done
                - integrate : the mesh will be used in order to compute some
                    integrals over the domain
                - FreeFEM : the mesh will be used in order to import/export mesh
                    file for FreeFEM

        """
        self.vertices = np.array(vertices)
        self.triangles = np.array(triangles)
        self.computations(domain, verbose =verbose, requirements=requirements)

    def add_convex_zone(self, zone_name: str, zone_vertices: list[PointType]) -> None:
        """Add a convex zone to the Mesh object.

        Args:
            zone_name (str) : the unique name designing the zone to add.
            zone_vertices (list[PointType]) :  the vertices of the outer boundary of
                the zone. The zone is supposed convex.

        """
        if(zone_name in self.zones):
            msg = f"Zone name {zone_name} already in use."
            raise NameError(msg)

        self.zones[zone_name] = {"boundary" : zone_vertices,
                                "triangles" : []}

    def compute_zones_triangles(self) -> None:
        """Compute the triangles included in each zone."""
        for index, center in enumerate(self.barycenters):
            for zone_name in self.zones:
                if self.in_zone(center, zone_name):
                    self.zones[zone_name]["triangles"].append(index)

    def in_zone(self, point: PointType, zone_name: str) -> bool:
        """Test if a given point is included in a given zone.

        Args:
            point (PointType) :  the point to test.
            zone_name (str) :  the name of the zone to test.

        """
        if zone_name not in self.zones:
            msg = f"The name {zone_name} does not correspond to a zone of the Mesh."
            raise NameError(msg)
        for i in range(1,len(self.zones[zone_name]["boundary"]) - 1):
            if(belong_triangle(point, [self.zones[zone_name]["boundary"][0],
                                       self.zones[zone_name]["boundary"][i],
                                       self.zones[zone_name]["boundary"][i+1]])):
                return True
        return False

    def export_mesh_msh(self, filename: str) -> None:
        """Export the data of the Mesh object in a .msh file.

        Args:
            filename (str) :  the complete path to the file.

        Raises:
            ImportError: requires meshio to be installed to work properly.

        """
        if not meshio:
            msg = "meshio must be installed in order to use .msh related methods."
            raise ImportError(msg)
        points = self.vertices
        cells = [
            ("line", np.array([self.edges[i] for i in np.concatenate([self.exit_edges,
                                                                      self.wall_edges])])),
            ("triangle", self.triangles),
        ]


        cell_data_dict = {
        "gmsh:geometrical" : [[4 for i in np.concatenate([self.exit_edges,
                                                          self.wall_edges])],
                              np.array([1 for t in self.triangles])],
        }

        mesh = meshio.Mesh(
            points,
            cells,
            point_data = {"gmsh:dim_tags" : np.array([[2,0] for P in self.vertices])},
            # Optionally provide extra data on points, cells, etc.
            cell_data = cell_data_dict,
            )
        mesh.write(
            filename,  # str, os.PathLike, or buffer/open file
            file_format="gmsh",  # optional if first argument is a path;
        )
        print("Saving mesh as ", filename)

    def export_mesh_msh_free_fem(self, filename: str) -> None:
        """Export the data of the Mesh object in a .msh file following the specific structure of the FreeFEM mesh files.

        Args:
            filename (str) :  the complete path to the file.

        """
        with Path(filename).open("w") as file:
            file.write(f"{len(self.vertices)}"
                        f" {len(self.triangles)}"
                        f" {len(self.exit_edges)+len(self.wall_edges)}\n")

            file.writelines([f"{point[0]} {point[1]} {self.vertex_flags[i]}\n"
                             for i,point in enumerate(self.vertices)])

            file.writelines([f"{triangle[0]+1} {triangle[1]+1} {triangle[2]+1} 0 \n"
                             for triangle in self.triangles])

            file.writelines([f"{self.edges[edge_index][0]+1}"
                             f" {self.edges[edge_index][1]+1} 98 \n"
                             for edge_index in self.exit_edges])

            file.writelines([f"{self.edges[edge_index][0]+1}"
                             f" {self.edges[edge_index][1]+1} 99 \n"
                             for edge_index in self.wall_edges])

        print("Saving mesh as ", filename)


    def generate_mesh_from_domain(self, domain: NonConvexDomain,
                                        dx: float,
                                        da: float=30,
                                        requirements : list[str] = DEFAULT_REQUIREMENTS,
                                        *,
                                        verbose :bool = False) -> None:
        """Compute a triangular mesh covering domain with the maximal length of an edge being set to dx.

        The computations are done using triangle (see https://www.cs.cmu.edu/~quake/triangle.html)

        Args:
            domain (NonConvexDomain): the domain object from which the mesh must be
                generated.
            dx (float): max area of a triangle (heavy computations when set low,
                computation time ~ 1/(dx)^2 ).
            da (float, optional): min angle inside a triangle (be careful, crash if set
                too high). Default = 30.
            verbose (bool, optional): displaying progression in console.
            requirements (list[str], optional): a list containing the computations that
                will be done using the mesh. The possible values are:

                - EikonalSolver : the mesh will be used in order to solve an eikonal
                    equation
                - LWRSolver : the mesh will be used in order to solve a scalar
                    conservation law
                - all : all the possible computations will be done
                - integrate : the mesh will be used in order to compute some
                    integrals over the domain
                - FreeFEM : the mesh will be used in order to import/export mesh
                    file for FreeFEM

        """
        domain_vertices = []
        domain_special_edges = []
        domain_vertices += domain.outer_vertices
        domain_special_edges += domain.outer_boundary

        for edge in domain.wall_edges:
            domain_special_edges.append([edge[0] +len(domain_vertices),
                                         edge[1] +len(domain_vertices)])
        domain_vertices += domain.wall_vertices

        if(len(domain.wall_holes_point)):
            domain_dict = {"vertices" : domain_vertices,
                           "segments" : domain_special_edges,
                           "holes" : domain.wall_holes_point,
                           }
        else:
            domain_dict = {"vertices" : domain_vertices,
                           "segments" : domain_special_edges,
                           }

        flag = "q"+str(da)+"pa"+str(dx)
        mesh_dict = tr.triangulate(domain_dict, flag)

        if "segments" not in mesh_dict:
            msg = "Corrupted domain, impossible to generate a mesh."
            raise ValueError(msg)
        if verbose :
            print(f"Mesh generated. Contains {len(mesh_dict['triangles'])} triangles")

        boundary_points = []
        for edge in mesh_dict["segments"]:
            for index in edge:
                if(index not in boundary_points):
                    boundary_points.append(index)

        self.boundary_points = np.array(boundary_points)

        self.vertices = mesh_dict["vertices"]
        self.triangles = mesh_dict["triangles"]
        self.computations(domain, verbose = verbose, requirements = requirements)

    def computations(self,  domain: NonConvexDomain = None,
                            requirements : list[str] = DEFAULT_REQUIREMENTS,
                            *,
                            verbose : bool = False) -> None:
        """Trigger all the computations required for the Mesh object.

        Args:
            domain (NonConvexDomain, optional): the domain object to use in order to
                recover the exits and the zones.
            verbose (bool, optional): displaying progression in console.
            requirements (list[str], optional): a list containing the computations that
                will be done using the mesh. The possible values are:

                - EikonalSolver : the mesh will be used in order to solve an eikonal
                    equation
                - LWRSolver : the mesh will be used in order to solve a scalar
                    conservation law
                - all : all the possible computations will be done
                - integrate : the mesh will be used in order to compute some
                    integrals over the domain
                - FreeFEM : the mesh will be used in order to import/export mesh
                    file for FreeFEM

        """
        if "all" in requirements:
            self.compute_edge_list()

            if domain is not None:
                self.set_exits_from_domain(domain)
                self.compute_vertex_flags(domain)
                for zone_name in domain.zones:
                    self.add_convex_zone(zone_name, domain.zones[zone_name])

            self.compute_pair_of_triangles_list()
            self.compute_outer_normals()
            self.min_cell_area = self.compute_cell_areas()
            self.max_edge_length = self.compute_edge_length()
            if(verbose):
                print("Minimal area for a triangle in the mesh : ", self.min_cell_area)
            self.compute_edge_length()
            self.compute_triangles_per_vertex()
            self.compute_zones_triangles()
        else:
            if("EikonalSolver" in requirements
               or "LWRSolver" in requirements
               or "FreeFEM" in requirements):
                self.compute_edge_list()

            if("EikonalSolver" in requirements or "LWRSolver" in requirements):
                self.set_exits_from_domain(domain)
                self.max_edge_length = self.compute_edge_length()
                for zone_name in domain.zones:
                    self.add_convex_zone(zone_name, domain.zones[zone_name])

            if("FreeFEM" in requirements):
                self.compute_vertex_flags(domain)

            if("LWRSolver" in requirements):
                self.compute_pair_of_triangles_list()
                self.compute_outer_normals()

            if("LWRSolver" in requirements or "integrate" in requirements):
                self.min_cell_area = self.compute_cell_areas()
                if(verbose):
                    print("Minimal area for a triangle in the mesh : ",
                          self.min_cell_area)

            if("EikonalSolver" in requirements or "integrate" in requirements):
                self.compute_triangles_per_vertex()

            if("integrate" in requirements):
                self.compute_zones_triangles()


    def compute_edge_list(self) -> None:
        """Fill two arrays concerning the edges of the mesh.

        The arrays are:

        - self.edges : cointaining all the edges of the mesh as a pair of vertex index.
        - self.triangles_with_edges : containing all the triangles of self.triangles,
            in the same order, represented as a triplet of edge index.

        Note:
            Required for:

                - LWRSolver
                - EikonalSolver
                - set_exits_from_domain

        """
        edges = []
        triangle_with_edge_list = []

        for triangle in self.triangles:
            edge1 = [triangle[0],triangle[1]]
            edge2 = [triangle[1],triangle[2]]
            edge3 = [triangle[2],triangle[0]]

            sides = [edge1,edge2,edge3]
            triangle_edge_coordinates = []

            for edge2 in sides:
                is_already_in = 0
                for n,edge in enumerate(edges):
                    if( (edge[0] == edge2[0] and edge[1] == edge2[1] )
                        or (edge[0] == edge2[1] and edge[1] == edge2[0]) ):
                        is_already_in = 1
                        triangle_edge_coordinates.append(n)
                        break
                if( not is_already_in):
                    edges.append(edge2)
                    triangle_edge_coordinates.append(len(edges)-1)

            triangle_with_edge_list.append(triangle_edge_coordinates)
        self.edges = np.array(edges)
        self.triangles_with_edges = np.array(triangle_with_edge_list)

    def set_exits_from_domain(self, domain: NonConvexDomain) -> None:
        """Configure the exits and the wall edges and vertices lists from the domain object passed in parameter.

        Args:
            domain (NonConvexDomain): the domain object containing information about
                the exits.

        Note:
            Required for:

                - LWRSolver
                - EikonalSolver

        """
        exit_vertices = {i for i in range(len(self.vertices))
                         for exit_e in domain.exit_list
                         if NonConvexDomain.belong_segment(self.vertices[i],exit_e)}
        self.exit_vertices = np.array(list(exit_vertices))

        exit_edges = []
        for i,edge in enumerate(self.edges):
            if(domain.has_exit_edge([self.vertices[edge[0]],self.vertices[edge[1]]])):
                exit_edges.append(i)
        self.exit_edges = np.array(exit_edges)

        if(len(exit_edges) == 0):
            msg = "Your mesh has no exit edge."
            raise ValueError(msg)

        wall_edges = []
        for index, edge in enumerate(self.edges):
            if( ( domain.has_wall_edge([self.vertices[edge[0]],
                                        self.vertices[edge[1]]])
                 and not domain.has_exit_edge([self.vertices[edge[0]],
                                               self.vertices[edge[1]]]) )
                or ( domain.has_outer_edge([self.vertices[edge[0]],
                                            self.vertices[edge[1]]])
                    and not domain.has_exit_edge([self.vertices[edge[0]],
                                                  self.vertices[edge[1]]]) )):
                wall_edges.append(index)
        self.wall_edges = np.array(wall_edges)

    def compute_vertex_flags(self, domain : NonConvexDomain,
                             flag_dict : dict = DEFAULT_MESH_DICT) -> None:
        """Compute the vertices flags necessary to export the Mesh as FreeFEM mesh file.

        Args:
            domain (NonConvexDomain, optional): the domain object to use in order to
                recover the exits and the zones.
            flag_dict (dict): dictionary prescribing the integer corresponding to the
                different types of vertex for an export as a FreeFEM .msh file.

        Note:
            Required for:

                - export_mesh_msh_free_fem

        """
        self.vertex_flags = []
        for point in self.vertices:
            if domain.has_exit_point(point):
                self.vertex_flags.append(flag_dict["exit"])
            elif domain.has_wall_point(point) or domain.has_boundary_point(point):
                self.vertex_flags.append(flag_dict["wall"])
            else:
                self.vertex_flags.append(flag_dict["domain"])


    def compute_triangles_per_vertex(self) -> None:
        """Compute the self.triangles_per_vertex list.

        Note:
            Required for:

                - EikonalSolver
                - integrateOverSquareBall

        """
        self.triangles_per_vertex : list = []

        for index in range(len(self.vertices)):
            triangle_list = []
            for triangle_index, triangle in enumerate(self.triangles):
                if index in triangle:
                    exclude = [int(P) for P in triangle if index != P]
                    triangle_list.append([triangle_index, exclude])
            self.triangles_per_vertex.append(triangle_list)


    def compute_pair_of_triangles_list(self) -> None:
        """Compute the lists pairs_of_triangles and boundary_edges_index.

        Note:
            Required for:

                - LWRSolver
                - checkGradientValidity

        """
        self.pairs_of_triangles=[]
        boundary_edges_index = []

        for index in range(len(self.edges)):
            pair_of_triangles = []
            number_of_loop = 0
            while(len(pair_of_triangles) < TWO
                  and number_of_loop < len(self.triangles_with_edges)):
                if(index in self.triangles_with_edges[number_of_loop]):
                    pair_of_triangles.append(number_of_loop)
                number_of_loop += 1
            self.pairs_of_triangles.append(pair_of_triangles)
            if(len(pair_of_triangles) < TWO):
                boundary_edges_index.append(index)
        self.boundary_edges_index = np.array(boundary_edges_index)

    def compute_outer_normals(self) -> None:
        """Compute the list outer_normal_vect_by_triangles.

        Note:
            Required for:

                - LWRSolver

        """
        outer_normal_vect_by_triangles = []

        def compute_outer_normal_unit_vect(a:PointType,b:PointType,c:PointType) -> PointType:
            norm = np.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2)
            x = (norm**2 + (c[0]-b[0])**2 + (c[1]-b[1])**2
                - (c[0]-a[0])**2 - (c[1]-a[1])**2 )/(2*norm)
            h = [ a[0] + (b[0]-a[0])*x/norm, a[1] + (b[1]-a[1])*x/norm ]
            sign = (c[0]-h[0])*(b[1]-a[1]) + (c[1]-h[1])*(a[0]-b[0])
            if(sign > 0):
                return [(a[1]-b[1])/norm,(b[0]-a[0])/norm]
            return [(b[1]-a[1])/norm,(a[0]-b[0])/norm]

        for index, triangle in enumerate(self.triangles_with_edges):
            normals = []
            normals = [compute_outer_normal_unit_vect(
                            self.vertices[self.edges[edgeindex][0]],
                            self.vertices[self.edges[edgeindex][1]],
                            self.vertices[point_index])
                       for edgeindex in triangle
                       for point_index in self.triangles[index]
                       if point_index not in self.edges[edgeindex]]
            if(len(normals) != 3):
                msg = (f"Corrupted mesh lists, found {len(normals)} sides"
                       " on one triangle.")
                raise ValueError(msg)
            outer_normal_vect_by_triangles.append(normals)
        self.outer_normal_vect_by_triangles = np.array(outer_normal_vect_by_triangles)

    def compute_cell_areas(self) -> float:
        """Compute the areas of each triangular cell.

        Raises an error if a triangle is degenerated.
        Computes also the barycenter of each triangle.

        Returns:
            float: the minimal area of a triangle in the mesh.

        Note:
            Required for:

                - integrate
                - LWRSolver

        """
        barycenters = []
        cell_areas = []
        min_area = 100
        for triangle in self.triangles:
            a = self.vertices[triangle[0]]
            b = self.vertices[triangle[1]]
            c = self.vertices[triangle[2]]

            barycenters.append([(a[0]+b[0]+c[0])/3,(a[1]+b[1]+c[1])/3])
            cell_areas.append(abs( (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]))/2)
            if(cell_areas[-1] < min_area):
                min_area = cell_areas[-1]
            elif(cell_areas[-1] == 0):
                msg = f"Degenerated mesh: the triangle [{a},{b},{c}] is of area 0."
                raise ValueError(msg)
        self.cell_areas = np.array(cell_areas)
        self.barycenters = np.array(barycenters)
        return(min_area)

    def compute_edge_length(self) -> float:
        """Compute the edge_length array.

        Note:
            Required for:

                - LWRSolver
                - EikonalSolver

        Returns:
            float: the maximal length of an edge in the mesh.

        """
        edge_length = []
        max_length = 0
        for edge in self.edges:
            a = self.vertices[edge[0]]
            b = self.vertices[edge[1]]

            new_len = np.sqrt( (b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]))
            edge_length.append(new_len)
            max_length = max( new_len, max_length)

        self.edge_length = np.array(edge_length)
        return max_length

    def get_limits(self) -> list[list[float]]:
        """Compute and return the extremal coordinates of the vertices of the mesh.

        Returns:
            list[list[float]]: the extremal coordinates of the vertices of the mesh as
                a list of two two-elements lists i.e. [[x_min,x_max],[y_min,y_max]].

        """
        x_min, x_max = self.vertices[0][0],self.vertices[0][0]
        y_min, y_max = self.vertices[0][1],self.vertices[0][1]
        for point in self.vertices:
            x_max = max(x_max,point[0])
            y_max = max(y_max,point[1])
            x_min = min(x_min,point[0])
            y_min = min(y_min,point[1])
        return [[x_min,x_max],[y_min,y_max]]

    def save_to_json(self, filename: str) -> None:
        """Save the mesh object in a .json file.

        Args:
            filename (str): the name of the file where the mesh will be saved.
                The extension .json is not needed in the filename.

        Raises:
            ImportError: the json module must be installed.

        """
        if not json:
            msg = "Module json not or wrongly installed. Needed for the json methods"
            raise ImportError(msg)

        mesh_dico = {"type": "triangular mesh"}
        mesh_dico["minCellArea"] = self.min_cell_area
        mesh_dico["maxEdgeLength"] = self.max_edge_length
        mesh_dico["vertices"] = self.vertices.tolist()
        mesh_dico["edges"] = self.edges.tolist()
        mesh_dico["triangles"] = self.triangles.tolist()

        mesh_dico["outerVertices"] = [self.vertices[i].tolist() for i in self.boundary_points]

        mesh_dico["exitVertices"] = self.exit_vertices.tolist()
        mesh_dico["exitEdges"] = self.exit_edges.tolist()
        mesh_dico["wallEdges"] = self.wall_edges.tolist()
        mesh_dico["boundaryPoints"] = self.boundary_points.tolist()
        mesh_dico["vertexFlags"] = self.vertex_flags
        mesh_dico["boundaryEdgesIndex"] = self.boundary_edges_index.tolist()

        mesh_dico["trianglesWithEdges"] = self.triangles_with_edges.tolist()
        mesh_dico["pairsOfTriangles"] = self.pairs_of_triangles
        mesh_dico["trianglesPerVertex"] = self.triangles_per_vertex

        mesh_dico["outerNormalVectByTriangles"] = self.outer_normal_vect_by_triangles.tolist()
        mesh_dico["cellAreas"] = self.cell_areas.tolist()
        mesh_dico["barycenters"] = self.barycenters.tolist()
        mesh_dico["edgeLength"] = self.edge_length.tolist()

        mesh_dico["zones"] = self.zones

        with Path(filename+"_mesh.json").open("w", encoding="utf-8") as f:
            json.dump(mesh_dico, f, ensure_ascii=False, indent=4)

    def load_from_json(self, filename: str) -> None:
        """Load the mesh object with the data from a .json file.

        Args:
            filename (str): the name of the file to load. The extension .json is needed
                in the filename.

        Raises:
            ImportError: the json module must be installed.

        """
        if not json:
            msg ="Module json not or wrongly installed. Needed for the json methods"
            raise ImportError(msg)

        with Path(filename).open("r") as f:
            data = json.load(f)
            self.min_cell_area : float = data["minCellArea"]
            self.max_edge_length : float = data["maxEdgeLength"]
            print("Minimal area for a triangle in the mesh : ", self.min_cell_area)

            self.vertices : ArrayLike = np.array(data["vertices"])
            self.edges : ArrayLike = np.array(data["edges"])
            self.triangles : ArrayLike = np.array(data["triangles"])

            self.exit_vertices : ArrayLike = np.array(data["exitVertices"])
            self.exit_edges : ArrayLike = np.array(data["exitEdges"])
            self.wall_edges : ArrayLike = np.array(data["wallEdges"])
            self.boundary_points : ArrayLike = np.array(data["boundaryPoints"])
            self.boundary_edges_index : ArrayLike = np.array(data["boundaryEdgesIndex"])
            self.vertex_flags : list = data["vertexFlags"]

            self.triangles_with_edges : ArrayLike = np.array(data["trianglesWithEdges"])
            self.pairs_of_triangles : list[list[int]] = data["pairsOfTriangles"]
            self.triangles_per_vertex : list[list] = data["trianglesPerVertex"]

            self.outer_normal_vect_by_triangles : ArrayLike = np.array(data["outerNormalVectByTriangles"])
            self.cell_areas : ArrayLike = np.array(data["cellAreas"])
            self.barycenters : ArrayLike = np.array(data["barycenters"])
            self.edge_length : ArrayLike = np.array(data["edgeLength"])

            self.zones :dict = data["zones"]
            print("Mesh successfully loaded.")



    def show(self, with_domain: NonConvexDomain=None, preference:str = "plotly") -> None:
        """Plot method for the Mesh object.

        Args:
            with_domain (NonConvexDomain, optional): domain object from which exits and
                walls can be recovered for plotting.
            preference (str, optional): set to "plotly" or "matplotlib" to chose the
                preferred plotting package. If only one package is installed the
                preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed
                and no preference is set, plotly is used.

        """
        if( go and (not plt or preference == "plotly") ):
            fig = go.Figure()
            if(with_domain):
                with_domain.add_plot(fig,ax=None,preference=preference)
            for triangle in self.triangles:
                fig.add_trace(go.Scatter(x=([self.vertices[i][0] for i in triangle]
                                            +[self.vertices[triangle[0]][0]]),
                                        y=([self.vertices[i][1] for i in triangle]
                                           +[self.vertices[triangle[0]][1]]),
                                fill="toself",
                                fillcolor="White",
                                mode="lines",
                                line={
                                    "color" : "Black",
                                    "width" : 1,
                                 }))
            fig.update_layout(yaxis={
                                        "scaleanchor" : "x",
                                        "scaleratio" : 1,
                                    })
            for edge in self.wall_edges:
                fig.add_shape(type="line",
                        x0=self.vertices[self.edges[edge][0]][0],
                        y0=self.vertices[self.edges[edge][0]][1],
                        x1=self.vertices[self.edges[edge][1]][0],
                        y1=self.vertices[self.edges[edge][1]][1],
                        line={
                            "color" : "LightSeaGreen",
                            "width" : 2,
                         })
            for edge in self.exit_edges:
                fig.add_shape(type="line",
                        x0=self.vertices[self.edges[edge][0]][0],
                        y0=self.vertices[self.edges[edge][0]][1],
                        x1=self.vertices[self.edges[edge][1]][0],
                        y1=self.vertices[self.edges[edge][1]][1],
                        line={
                            "color" : "Red",
                            "width" : 2,
                         })

            fig.show()
        elif(plt):
            fig, ax = plt.subplots()
            self.add_plot(fig,ax,preference)
            for edge in self.wall_edges:
                path = [    (self.vertices[self.edges[edge][0]][0],
                             self.vertices[self.edges[edge][0]][1]),
                            (self.vertices[self.edges[edge][1]][0],
                             self.vertices[self.edges[edge][1]][1]) ]
                ax.add_patch(patches.Polygon(path, edgecolor="black",linewidth=1))

            for exit_edge in self.exit_edges:
                path = [    (self.vertices[self.edges[exit_edge][0]][0],
                             self.vertices[self.edges[exit_edge][0]][1]),
                            (self.vertices[self.edges[exit_edge][1]][0],
                             self.vertices[self.edges[exit_edge][1]][1]) ]
                ax.add_patch(patches.Polygon(path, edgecolor="red",linewidth=2))
            limits = self.get_limits()
            ax.set_xlim(limits[0][0],limits[0][1])
            ax.set_ylim(limits[1][0],limits[1][1])
            plt.axis("equal")
            plt.show()
        else:
            msg = ("No plotting module found. Try installing plotly or matplotlib"
                   " if you want to use show methods")
            raise ImportError(msg)

    def add_plot(self, fig: FigureType, ax: AxesType = None, preference: str = "plotly") -> None:
        """Non-blocking plot method for the mesh object.

        Args:
            fig (matplotlib.Figure or plotly.graph_objects.Figure): the instance of
                Figure to which the plot must be added. Depends on the library used.
            ax (matplotlib.Axes or None): the instance of Axes to use if matplotlib is
                used.
            preference (str): set to "plotly" or "matplotlib" to chose the preferred
                plotting package. If only one package is installed the preference is
                ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed for this method to
                work. If both are installed and no preference is set, plotly is used.

        """
        if( go and (not plt or preference == "plotly") ):
            for triangle in self.triangles:
                fig.add_trace(go.Scatter(x=([self.vertices[i][0] for i in triangle]
                                            +[self.vertices[triangle[0]][0]]),
                                        y=([self.vertices[i][1] for i in triangle]
                                           +[self.vertices[triangle[0]][1]]),
                                fill="toself",
                                fillcolor="White",
                                mode="lines",
                                line={
                                    "color" : "Black",
                                    "width" : 0.5,
                                 }))
        elif(plt):
            for triangle in self.triangles:
                triangle_fig = patches.Polygon([(self.vertices[i][0],
                                                 self.vertices[i][1])
                                                for i in triangle],
                                               edgecolor="black",
                                               facecolor="white")
                ax.add_patch(triangle_fig)

        else:
            msg = ("No plotting module found. Try installing plotly or matplotlib"
                   " if you want to use show methods")
            raise ImportError(msg)


class CellValueMap:
    """An object to represent a cell valued map that corresponds to a function which is constant on the triangles of a mesh (mostly densities).

    Attributes
    ----------

    Attributes:
        Mesh (Mesh): the corresponding Mesh object on which the map is defined.
        values (list[float]): the values on each triangle as a list ordered in the same
            way as Mesh.triangles.

    Methods
    -------

    """

    def __init__(self, mesh: Mesh) -> None:
        self.mesh : Mesh = mesh
        self.values : list[float] = [0 for _ in self.mesh.triangles]

    def generate_random(self,variability :float = 0.4, mean : float= 0.23) -> None:
        """Generate random values for the map on each triangle of the mesh.

        The random value follows a uniform law on
        [mean -variability/2, mean + variability/2].

        Args:
            variability (float): the length of the interval of random values.
            mean (float): the center of the interval of random values.

        """
        self.values = [ mean + variability*(alea.random()-0.5)
                       for _ in self.mesh.triangles]

    def set_constant(self, value: float) -> None:
        """Set the value to "value" on every triangle.

        Args:
            value (float): the value to set on each triangle.

        """
        for i in range(len(self.values)):
            self.values[i] = value

    def __len__(self) -> int:
        return len(self.values)

    def __getitem__(self, index:int) -> float:
        return self.values[index]

    def __setitem__(self, index:int, value:float) -> None:
        self.values[index] = value

    def __add__(self, other: Self|float ) -> Self:
        if isinstance(other, CellValueMap):
            result = CellValueMap(self.mesh)
            for index in range(len(self.mesh.triangles)):
                result[index] = self[index] + other[index]
            return result
        if isinstance(other, (float, int)):
            result = CellValueMap(self.mesh)
            for index in range(len(self.mesh.triangles)):
                result[index] = self[index] + other
            return result

        msg = ("A CellValueMap can be summed only with another CellValueMap"
               " or a float or an integer.")
        raise TypeError(msg)

    def __mul__(self, other:Self|float) -> Self:
        if isinstance(other, CellValueMap):
            result = CellValueMap(self.mesh)
            for index in range(len(self.mesh.triangles)):
                result[index] = self[index]*other[index]
            return result
        if isinstance(other, (float, int)):
            result = CellValueMap(self.mesh)
            for index in range(len(self.mesh.triangles)):
                result[index] = self[index]*other
            return result

        msg = ("A CellValueMap can only be multiplied with another CellValueMap"
                " or a float or an integer.")
        raise TypeError(msg)

    def integrate(self) -> float:
        """Compute the integral of the map over the whole domain.

        Returns:
            float: the integral of the CellValueMap over the whole domain.

        """
        return sum([self.values[i]*self.mesh.cell_areas[i]
                    for i in range(len(self.mesh.triangles))])

    def set_constant_circle(self, center:list[PointType], radius:float, value:float) -> None:
        """Set to the given value all the cells for which the barycenter is inside the disk of given center and radius.

        Args:
            center (list[PointType]): the center of the disk
            radius (float): the radius of the disk
            value (float): the value to set on the corresponding cells.

        """
        for index in range(len(self.mesh.triangles)):
            if(((self.mesh.barycenters[index][0] - center[0])**2
                + (self.mesh.barycenters[index][1] - center[1])**2 )
               <= radius**2):
                self.values[index] = value

    def convolution_over_square_ball(self, radius:float, conv_func:Callable[[float], float]) -> list:
        r"""Compute the convolution of the CellValueMap with the conv_func function on the support defined by the square of given radius (infinity-norm ball).

        Args:
            radius (float): the radius :math:`r` of the convolution support.
            conv_func (function, float -> float): the convolution function in the form
                :math:`F(\rho(y),|x_1-y_1|,|x_2-y_2|)` where :math:`\rho` is the
                cellValueMap and :math:`x = (x_1,x_2)` is the vertex where the
                convolution is computed. The computed quantity is then,
                for any :math:`x = (x_1,x_2)`,

                .. math::
                    \iint_{[-r/2, r/2]^2} F(\rho(x+y), |y_1|,|y_2|) d y_1 d y_2.

        Returns:
            list: a list containing the computed value of the convolution for each
            vertex of the mesh ordered in the same way as the vertices of the
            Mesh object.

        """
        def recursive_integral(center:list[float], rad:float, visited:list[int], index:int) -> float:
            visited.append(index)
            dist_x = np.abs(self.mesh.barycenters[index][0] - center[0])
            dist_y = np.abs(self.mesh.barycenters[index][1] - center[1])
            if(dist_x > rad or dist_y > rad):
                return 0.0
            conv_sum = conv_func(self.values[index],dist_x,dist_y)*self.mesh.cell_areas[index]
            for vertex in self.mesh.triangles[index]:
                for neighbor_triangle in self.mesh.triangles_per_vertex[vertex]:
                    if neighbor_triangle[0] not in visited:
                        conv_sum += recursive_integral(center,rad, visited, neighbor_triangle[0])
            return conv_sum

        return [recursive_integral(self.mesh.vertices[vertexIndex],radius,[],
                                   self.mesh.triangles_per_vertex[vertexIndex][0][0])
                for vertexIndex in range(len(self.mesh.vertices))]

    def fit_averaged_map(self, other: Self) -> None:
        """Fit a cellValueMap over another cellValueMap when the Mesh are different but the domains are the same by averaging over each triangle.

        Args:
            other (CellValueMap): the other map to fit onto the self map.

        """
        for i in range(len(self.values)):
            new_len_cell = self.mesh.points[i+1] - self.mesh.points[i]
            start = 0
            while(other.Mesh.points[start] < self.mesh.points[i]):
                start += 1
            if(other.Mesh.points[start] > self.mesh.points[i]):
                start -= 1
            end = start+1
            while(other.Mesh.points[end] < self.mesh.points[i+1]):
                end += 1
            j = start
            average = 0
            while(j <= end -1):
                average += (other.values[j]
                            *(min(other.Mesh.points[j+1], self.mesh.points[i+1])
                              - max(other.Mesh.points[j], self.mesh.points[i])))
                j += 1
            self.values[i] = average/new_len_cell

    def show(self, preference:str = "plotly") -> None:
        """Plot method for the CellValueMap object.

        Args:
            preference (str, optional): set to "plotly" or "matplotlib" to chose the
                preferred plotting package. If only one package is installed the
                preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed
                and no preference is set, plotly is used.


        """
        if( go and (not plt or preference == "plotly") ):
            fig = go.Figure()
            for j,triangle in enumerate(self.mesh.triangles):
                fig.add_trace(go.Scatter(x=([self.mesh.vertices[i][0] for i in triangle]
                                            +[self.mesh.vertices[triangle[0]][0]]),
                                        y=([self.mesh.vertices[i][1] for i in triangle]
                                           +[self.mesh.vertices[triangle[0]][1]]),
                                fill="toself",
                                hoverinfo = "none",
                                showlegend = False,
                                mode="none",
                                fillcolor =("rgb("
                                            +str(int(255*min(1,max(self.values[j],0))))
                                            +",0,0)"),
                                ))
            fig.update_layout(yaxis={
                                        "scaleanchor" : "x",
                                        "scaleratio" : 1,
                                    })
            fig.show()
        elif(plt):
            fig, ax = plt.subplots()

            col = collections.PolyCollection([[ (self.mesh.vertices[i][0],
                                                 self.mesh.vertices[i][1])
                                               for i in triangle]
                                              for triangle in self.mesh.triangles])
            col.set_cmap(cm.viridis)
            col.set_clim([0, 1])
            rgcol = col.set_array(self.values)

            ax.add_collection(col)
            fig.colorbar(rgcol, ax=ax, label="density")

            limits = self.mesh.get_limits()
            ax.set_xlim(limits[0][0],limits[0][1])
            ax.set_ylim(limits[1][0],limits[1][1])
            plt.axis("equal")
            plt.show()
        else:
            msg = ("No plotting module found. Try installing plotly or matplotlib"
                   " if you want to use show methods")
            raise ImportError(msg)

    def get_scatter(self, preference: str = "plotly") -> list[FigureType]:
        """Non-blocking plot method for the cellValueMap object.

        Args:
            preference (str, optional): set to "plotly" or "matplotlib" to chose the
                preferred plotting package. If only one package is installed the
                preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed
                and no preference is set, plotly is used.

        """
        if( go and (not plt or preference == "plotly") ):
            list_plots = []
            for j,triangle in enumerate(self.mesh.triangles):
                list_plots.append(
                            go.Scatter(x=([self.mesh.vertices[i][0] for i in triangle]
                                       +[self.mesh.vertices[triangle[0]][0]]),
                                        y=([self.mesh.vertices[i][1] for i in triangle]
                                           +[self.mesh.vertices[triangle[0]][1]]),
                                fill="toself",
                                hoverinfo = "none",
                                showlegend = False,
                                mode="none",
                                fillcolor =("rgb("
                                            +str(int(255*min(1,max(self.values[j],0))))
                                            +",0,0)"),
                                ))
            return list_plots
        if(plt):
            print("Useless method for matplotlib.")
            return []

        msg = ("No plotting module found. Try installing plotly or matplotlib"
               " if you want to use show methods")
        raise ImportError(msg)


class VertexValueMap:
    """An object to represent a vertex valued map that corresponds to a function which is affine on the triangles, defined by its values on the vertices (often a potential).

    Attributes
    ----------

    Attributes:
        Mesh (Mesh): the corresponding Mesh object on which the map is defined.
        values (list[float]): the values on each vertex as a list ordered in the same
            way as Mesh.vertices.

    Methods
    -------

    """

    def __init__(self, mesh:Mesh) -> None:
        self.mesh = mesh
        self.values = np.array([0 for _ in self.mesh.vertices])

    def generate_random(self, variability:float = 0.5) -> None:
        """Generate random values for the maps distributed as 0.5 + variability*X where X is uniform on [-0.5,0.5].

        Args:
            variability (float, optional): the range of the uniform distribution.

        """
        self.values = np.array([ 0.5 + variability*(alea.random()-0.5)
                                for _ in self.mesh.vertices])

    def __len__(self) -> int:
        return len(self.values)

    def __getitem__(self, index:int) -> float:
        return self.values[index]

    def __setitem__(self, index:int, value:float) -> None:
        self.values[index] = value

    def __add__(self,other:Self) -> Self:
        new_one = VertexValueMap(self.mesh)
        new_one.values = self.values + other.values
        return new_one

    def __str__(self) -> str:
        return(str(self.values))

    def set_infinity(self) -> None:
        """Set all values of the map to *float("inf")*."""
        self.values = [float("inf") for _ in self.mesh.vertices]


    def show_3d(self, preference: str = "plotly") -> None:
        """Display a 3D plot of the vertex valued map.

        Note:
            Only available with plotly at the moment.

        Args:
            preference (str, optional): set to "plotly" or "matplotlib" to chose the
                preferred plotting package. If only one package is installed the
                preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed
                and no preference is set, plotly is used.

        """
        if( go and (not plt or preference == "plotly") ):
            fig = go.Figure()
            fig.add_trace(go.Mesh3d(x=[P[0] for P in self.mesh.vertices],
                                    y = [P[1] for P in self.mesh.vertices],
                                    z = self.values,
                                    opacity=1,
                                    color="rgba(244,22,100,0.6)"))
            fig.show()
        else:
            msg = ("No plotting module found. Try installing plotly or matplotlib"
                   " if you want to use show methods")
            raise ImportError(msg)

    def add_3d_plot(self, fig:FigureType, color:list[float]=DEFAULT_COLOR_RGBA, preference:str ="plotly") -> None:
        """Non-blocking method that adds the 3D plot to a given `Figure` object.

        Note:
            Only available with plotly at the moment.

        Args:
            fig (matplotlib.pyplot.Figure or plotly.graph_objects.Figure): the instance
                of Figure to which the plot must be added. Depends on the library used.
            color (list[float]): color of the plot, the RGBA code.
            preference (str, optional): set to "plotly" or "matplotlib" to chose the
                preferred plotting package. If only one package is installed the
                preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed
                and no preference is set, plotly is used.

        """
        if( go and (not plt or preference == "plotly") ):
            fig.add_trace(go.Mesh3d(x=[P[0] for P in self.mesh.vertices],
                                    y = [P[1] for P in self.mesh.vertices],
                                    z = self.values,
                                    opacity=1,
                                    color="rgba("+str(color[0])+","+str(color[1])+","+str(color[2])+","+str(color[3])+")"))
        else:
            msg = ("No plotting module found. Try installing plotly or matplotlib"
                   " if you want to use show methods")
            raise ImportError(msg)


    def show(self, colorscale_name: str = "viridis", preference: str = "plotly", *, grid:bool = False) -> None:
        """Display the vertex valued map as a colorscale scatter plot in 2D.

        Note:
            Only available with plotly at the moment.

        Args:
            grid (bool, keyword-only): determines if the mesh is plotted or not.
            colorscale_name (str, optional): name of the colorscale to use for the plot.
            preference (str, optional): set to "plotly" or "matplotlib" to chose the
                preferred plotting package. If only one package is installed the
                preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed
                and no preference is set, plotly is used.

        """
        if( go and (not plt or preference == "plotly") ):
            fig = go.Figure()
            if(grid):
                self.mesh.add_plot(fig,ax=None,preference=preference)
            fig.add_trace(go.Scatter(x=[P[0] for P in self.mesh.vertices],
                                    y = [P[1] for P in self.mesh.vertices],
                            hoverinfo = "none",
                            showlegend = False,
                            mode="markers",
                            marker = {
                                    "color" : self.values,
                                    "colorscale" : colorscale_name,
                                }))
            fig.update_layout(yaxis={
                                        "scaleanchor" : "x",
                                        "scaleratio" : 1,
                                    })
            fig.show()
        else:
            msg = ("No plotting module found. Try installing plotly or matplotlib"
                   " if you want to use show methods")
            raise ImportError(msg)


    def compute_gradient_flow(self, normalization: Callable[[float,float],float] = (lambda x,y : np.sqrt(x**2 + y**2)), *, normalize:bool = True) -> ArrayLike:
        """Compute and return the gradients of the function defined by the vertex valued map. The gradient is computed on each triangle.

        Args:
            normalize (bool, keyword-only): determines if the gradients should be
                renormalized.
            normalization (function, (float, float) -> float): the norm to use for the
                renormalization if `normalize=True`.

        Raises:
            ValueError: if at least one of the triangles is degenerated.

        Returns:
            ArrayLike: an array containing all the gradients ordered as the list
                Mesh.triangles.

        """
        list_triangles_grads = []

        for triangle in self.mesh.triangles:
                det =  ((self.mesh.vertices[triangle[1]][0]
                         - self.mesh.vertices[triangle[0]][0])
                        *(self.mesh.vertices[triangle[2]][1]
                          - self.mesh.vertices[triangle[0]][1])
                        - (self.mesh.vertices[triangle[2]][0]
                           - self.mesh.vertices[triangle[0]][0])
                        *(self.mesh.vertices[triangle[1]][1]
                          - self.mesh.vertices[triangle[0]][1]))
                if(det == 0):
                    msg = (f"The triangle [{self.mesh.vertices[triangle[0]]},"
                            f"{self.mesh.vertices[triangle[1]]},"
                            f"{self.mesh.vertices[triangle[2]]}] is degenerated.")
                    raise ValueError(msg)
                vec_x = ( (self.values[triangle[0]] - self.values[triangle[2]])
                         *(self.mesh.vertices[triangle[1]][1]
                           - self.mesh.vertices[triangle[0]][1])
                        + (self.values[triangle[1]] - self.values[triangle[0]])
                        *(self.mesh.vertices[triangle[2]][1]
                          - self.mesh.vertices[triangle[0]][1]) )/det
                vec_y = ( (self.values[triangle[0]] - self.values[triangle[1]])
                         *(self.mesh.vertices[triangle[2]][0]
                           - self.mesh.vertices[triangle[0]][0])
                        + (self.values[triangle[2]] - self.values[triangle[0]])
                        *(self.mesh.vertices[triangle[1]][0]
                          - self.mesh.vertices[triangle[0]][0]) )/det
                list_triangles_grads.append([-vec_x/normalization(vec_x, vec_y),
                                             -vec_y/normalization(vec_x, vec_y)]
                                            if normalize else [-vec_x,-vec_y])
        return np.array(list_triangles_grads)

    def compute_vertex_gradient_flow(self, normalization: Callable[[float,float],float] = (lambda x,y : np.sqrt(x**2 + y**2)), *,  normalize:bool = True) -> ArrayLike:
        """Compute and return the gradients of the function defined by the vertex valued map.

        The gradient is computed at each vertex as a weighted mean of the differences
        with the neighbouring vertices.

        Note:
            This method is less precise than `compute_gradient_flow` in general.
            However due to its non-local construction, the gradient flow obtained with
            this method tends to be more regular.

        Args:
            normalize (bool, keyword-only): determines if the gradients should be
                renormalized.
            normalization (function, float, float -> float): the norm to use for the
                renormalization if `normalize=True`.

        Returns:
            ArrayLike: an array containing all the gradients ordered as the
                list Mesh.triangles.

        """
        vertex_triangles_grads = []

        for vertex in range(len(self.mesh.vertices)):
            mean_vect = np.array([0.0,0.0])
            treated_vertices = []
            for i in range(len(self.mesh.triangles_per_vertex[vertex])):
                for point in self.mesh.triangles_per_vertex[vertex][i][1]:
                    if point not in treated_vertices:
                        treated_vertices.append(point)
                        mean_vect += (((self.values[point] - self.values[vertex])
                                      /np.sqrt((self.mesh.vertices[point][0]
                                                -self.mesh.vertices[vertex][0])**2
                                               + (self.mesh.vertices[point][1]
                                                  -self.mesh.vertices[vertex][1])**2))
                                      *(self.mesh.vertices[point]-self.mesh.vertices[vertex]))
            vertex_triangles_grads.append([
                        mean_vect[0]/normalization(mean_vect[0], mean_vect[1]),
                        mean_vect[1]/normalization(mean_vect[0], mean_vect[1])]
                                  if normalize else mean_vect/len(treated_vertices))
        mean_triangle_grads = []
        for triangle in self.mesh.triangles:
            vec_x = sum([vertex_triangles_grads[i][0] for i in triangle])
            vec_y = sum([vertex_triangles_grads[i][1] for i in triangle])
            mean_triangle_grads.append([-vec_x/normalization(vec_x, vec_y),
                                        -vec_y/normalization(vec_x, vec_y)]
                                       if normalize else [-vec_x,-vec_y])
        return np.array(mean_triangle_grads)


    def show_vector_field(self, normalization: Callable[[float,float],float] = (lambda x,y : np.sqrt(x**2 + y**2)), preference: str = "plotly", *, normalize:bool = True) -> None:
        """Display the vector field corresponding to the gradient flow (obtained with `compute_gradient_flow`) of the vertex valued map.

        Note:
            Only available with plotly at the moment.

        Args:
            normalize (bool, optional): determines if the gradients should be
                renormalized.
            normalization (function, (float, float) -> float): the norm to use for the
                renormalization if `normalize=True`.
            preference (str, optional): set to "plotly" or "matplotlib" to chose the
                preferred plotting package. If only one package is installed the
                    preference is ignored.

        Raises:
            ImportError: plotly or matplotlib should be installed. If both are installed
                and no preference is set, plotly is used.

        """
        if( go and (not plt or preference == "plotly") ):
            vector_field = self.compute_gradient_flow(normalize = True,
                                                      normalization = normalization)
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=[P[0] for P in self.mesh.vertices],
                                    y = [P[1] for P in self.mesh.vertices],
                            hoverinfo = "none",
                            showlegend = False,
                            mode="markers",
                            marker = {
                                    "color" : self.values,
                                    "colorscale" : "viridis",
                                }))
            fig_quiv = ff.create_quiver([(self.mesh.vertices[triangle[0]][0]
                                          +self.mesh.vertices[triangle[1]][0]
                                          +self.mesh.vertices[triangle[2]][0])/3
                                         for triangle in self.mesh.triangles],
                                        [(self.mesh.vertices[triangle[0]][1]
                                          +self.mesh.vertices[triangle[1]][1]
                                          +self.mesh.vertices[triangle[2]][1])/3
                                         for triangle in self.mesh.triangles],
                                        [V[0] for V in vector_field],
                                        [V[1] for V in vector_field])
            fig.add_traces(fig_quiv.data)
            fig.update_layout(yaxis={
                                        "scaleanchor" : "x",
                                        "scaleratio" : 1,
                                    })
            fig.show()
        else:
            msg = ("No plotting module found. Try installing plotly or matplotlib"
                   " if you want to use show methods")
            raise ImportError(msg)
