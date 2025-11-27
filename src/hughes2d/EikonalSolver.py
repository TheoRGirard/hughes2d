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

from collections.abc import Callable

import numpy as np

from hughes2d.Mesh2D import CellValueMap, Mesh, VertexValueMap

EMPTY_MAP = []

class EikoSolver:
    r"""An object wrapping together all the methods used to solve the eikonal equation in a discontinuous context.

    We use the following notations for the eikonal equation:

    .. math::
        |\\nabla u | = c(\\rho).


    Attributes
    ----------

    Attributes:
        mesh (Mesh): the mesh on which the solution to the eikonal equation will be
            approximated.
        density (CellValueMap): a function constant on the triangle of the mesh
            representing the density of pedestrians.
        fieldValues (VertexValueMap): the object corresponding to the approximated
            solution of the eikonal equation.
        cost (function, float -> float): the function *c* that must be applied to the
            density map. If we set :math:`c: \\rho \\mapsto 1` the solution of the
            eikonal equation will correspond to the length of the shortest path
            to the exits, without taking the density into account.
        opt (dict): the dictionary containing all the parameters of computation:

            ================= ====== ================ ===================================================================================================================
            key               type   possible values  description
            ================= ====== ================ ===================================================================================================================
            "method"           str    "FME" or "FMT"   corresponds to the computational method to use
            "NarrowBandDepth"  int    1 or 2           the maximal degree of the neighbours explored at each iteration (only relevant if method="FMT")
            "constrained"      bool   True or False    determines if the gradient computed must be constrained inside the triangle or not (only relevant if method="FMT")
            "debugging"        bool   True or False    debugging option, for verbose outputs and detailled prints
            ================= ====== ================ ===================================================================================================================

    Note:
        The best combination of options with respect to the quality of the approximation
        appears to be::

            {
                "method" : "FMT",
                "NarrowBandDepth" : 2,
                "constrained" : True
            }

        For further details about the algorithms and their comparison, we defer to https://hal.science/tel-05275359v1

    Methods
    -------
    """

    def __init__(self, mesh: Mesh, density_map: CellValueMap = EMPTY_MAP, cost_function: Callable[[float], float] = (lambda x: 1+x), opt:dict = {}) -> None:
        self.mesh = mesh
        if(density_map != EMPTY_MAP):
            self.density = density_map
        else:
            self.density = CellValueMap(self.mesh)

        self.field_values = VertexValueMap(self.mesh)

        self.opt = opt
        if("method" not in self.opt):
            self.opt["method"] = "FMT"

        if("NarrowBandDepth" not in self.opt):
            self.opt["NarrowBandDepth"] = 2
        if("constrained" not in self.opt):
            self.opt["constrained"] = True

        if("debugging" not in self.opt):
            self.opt["debugging"] = False

        self.cost = cost_function

    def update_density(self, density: CellValueMap) -> None:
        """Update the density of the solver with the new density passed as a parameter.

        Args:
            density (CellValueMap): the new density to update in the solver.

        """
        self.density = density

    def compute_field_unconstrained_dep_2(self) -> None:
        """Compute the approximated solution of the eikonal equation and stores the approximation in "fieldValues".

        The method used is the FMT algorithm with no constraint on the gradient and
        vertices at distance one or two are considered neighbour.
        i.e. constrained = False and NarrowBandDepth = 2 see "opt"
        """
        state_map = [0 for i in range(len(self.mesh.vertices))]

        visited_dict = {}
        visited = []
        validated = []
        self.field_values.set_infinity()

        #Dirichlet conditions on the exits
        for index in self.mesh.exit_vertices:
            self.field_values[index] = 0
            validated.append(index)
            state_map[index] = 1

        #Construction of the narrow band + computation of the potential values
        for index in validated:
            for [triangleindex, triangle_without_index] in self.mesh.triangles_per_vertex[index]:

                if(state_map[triangle_without_index[1]] == 0):
                    state_map[triangle_without_index[1]] = 2

                if(state_map[triangle_without_index[0]] == 0):
                    state_map[triangle_without_index[0]] = 2

                if(state_map[triangle_without_index[0]] == 1
                   and state_map[triangle_without_index[1]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_unconstrained(
                                    self.mesh.vertices[triangle_without_index[1]],
                                    self.mesh.vertices[triangle_without_index[0]],
                                    self.mesh.vertices[index],
                                    self.field_values[triangle_without_index[0]],
                                    self.field_values[index],
                                    self.cost(self.density[triangleindex]))

                    if(triangle_without_index[1] in visited_dict):
                        if((triangleindex in visited_dict[triangle_without_index[1]])
                            and (visited_dict[triangle_without_index[1]][triangleindex]
                                 == self.field_values[triangle_without_index[1]])):
                            self.field_values[triangle_without_index[1]] = min([
                                   visited_dict[triangle_without_index[1]][tridex]
                                   for tridex in visited_dict[triangle_without_index[1]]
                                   if tridex != triangleindex]+[potential_value])
                            visited = self.add_in_order_by_field_value(visited,
                                                                  triangle_without_index[1])
                        visited_dict[triangle_without_index[1]][triangleindex] = potential_value
                    else:
                        visited_dict[triangle_without_index[1]] = {}
                        visited_dict[triangle_without_index[1]][triangleindex] = potential_value


                    if potential_value < self.field_values[triangle_without_index[1]]:
                        self.field_values[triangle_without_index[1]] = potential_value
                        visited = self.add_in_order_by_field_value(visited,
                                                              triangle_without_index[1])


                elif(state_map[triangle_without_index[1]] == 1
                     and state_map[triangle_without_index[0]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_unconstrained(
                                                    self.mesh.vertices[triangle_without_index[0]],
                                                    self.mesh.vertices[triangle_without_index[1]],
                                                    self.mesh.vertices[index],
                                                    self.field_values[triangle_without_index[1]],
                                                    self.field_values[index],
                                                    self.cost(self.density[triangleindex]))
                    if(triangle_without_index[0] in visited_dict):
                        if((triangleindex in visited_dict[triangle_without_index[0]])
                            and (visited_dict[triangle_without_index[0]][triangleindex]
                                 == self.field_values[triangle_without_index[0]])):
                            self.field_values[triangle_without_index[0]] = min([
                                 visited_dict[triangle_without_index[0]][tridex]
                                 for tridex in visited_dict[triangle_without_index[0]]
                                 if tridex != triangleindex]+[potential_value])
                            visited = self.add_in_order_by_field_value(visited,
                                                                  triangle_without_index[0])
                        visited_dict[triangle_without_index[0]][triangleindex] = potential_value
                    else:
                        visited_dict[triangle_without_index[0]] = {}
                        visited_dict[triangle_without_index[0]][triangleindex] = potential_value

                    if potential_value < self.field_values[triangle_without_index[0]]:
                        self.field_values[triangle_without_index[0]] = potential_value
                        visited = self.add_in_order_by_field_value(visited,
                                                              triangle_without_index[0])

                elif(state_map[triangle_without_index[1]] == 2
                     and state_map[triangle_without_index[0]] == 2):

                    potential_value = (self.field_values[index]
                                        + self.cost(self.density[triangleindex])
                                            *EikoSolver.compute_height_length(self.mesh.vertices[index],
                                                                            self.mesh.vertices[triangle_without_index[0]],
                                                                            self.mesh.vertices[triangle_without_index[1]]))

                    for i in range(2):
                        if( triangle_without_index[i] not in visited_dict):
                            visited_dict[triangle_without_index[i]] = {}

                        visited_dict[triangle_without_index[i]][triangleindex] = potential_value
                        if potential_value < self.field_values[triangle_without_index[i]]:
                            self.field_values[triangle_without_index[i]] = potential_value
                            visited = self.add_in_order_by_field_value(visited,
                                                                  triangle_without_index[i])


        #Main loop : validation of a vertex and recomputation of the narrow band
        last_validated = validated[0]
        num_vertices = len(state_map)

        while(len(validated) < num_vertices):

            #Validation of a vertex
            if(len(visited) == 0):
                if(self.opt["debugging"]):
                    print("Completion : ", len(validated), "/", num_vertices)
                msg = "The fast marching ended too fast. Is the Mesh connected ?"
                raise ValueError(msg)
            if(self.field_values[visited[0]] == float("inf")):
                if(self.opt["debugging"]):
                    print("Infinity value detected at vertex (",
                          self.mesh.vertices[visited[0]][0],",",
                          self.mesh.vertices[visited[0]][1], ")" )
                msg = "Gradient computation impossible."
                raise ValueError(msg)

            validated.append(visited[0])
            last_validated = visited[0]
            visited = visited[1:]
            state_map[last_validated] = 1

            #recomputation of the narrow band
            for [triangleindex, triangle_without_index] in self.mesh.triangles_per_vertex[last_validated]:

                if(state_map[triangle_without_index[1]] == 0):
                    state_map[triangle_without_index[1]] = 2

                if(state_map[triangle_without_index[0]] == 0):
                    state_map[triangle_without_index[0]] = 2

                if(state_map[triangle_without_index[0]] == 1
                   and state_map[triangle_without_index[1]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_unconstrained(
                                            self.mesh.vertices[triangle_without_index[1]],
                                            self.mesh.vertices[triangle_without_index[0]],
                                            self.mesh.vertices[last_validated],
                                            self.field_values[triangle_without_index[0]],
                                            self.field_values[last_validated],
                                            self.cost(self.density[triangleindex]))
                    if(triangle_without_index[1] in visited_dict): #Update visited dict
                        if((triangleindex in visited_dict[triangle_without_index[1]])
                            and(visited_dict[triangle_without_index[1]][triangleindex]
                                == self.field_values[triangle_without_index[1]])):
                            self.field_values[triangle_without_index[1]] = min([
                                visited_dict[triangle_without_index[1]][tridex]
                                for tridex in visited_dict[triangle_without_index[1]]
                                if tridex != triangleindex]+[potential_value])
                            visited = self.add_in_order_by_field_value(visited,
                                                                  triangle_without_index[1])
                        visited_dict[triangle_without_index[1]][triangleindex] = potential_value
                    else:
                        visited_dict[triangle_without_index[1]] = {}
                        visited_dict[triangle_without_index[1]][triangleindex] = potential_value

                    if potential_value < self.field_values[triangle_without_index[1]]:
                        self.field_values[triangle_without_index[1]] = potential_value
                        visited = self.add_in_order_by_field_value(visited,
                                                              triangle_without_index[1])


                elif(state_map[triangle_without_index[1]] == 1
                     and state_map[triangle_without_index[0]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_unconstrained(
                                        self.mesh.vertices[triangle_without_index[0]],
                                        self.mesh.vertices[triangle_without_index[1]],
                                        self.mesh.vertices[last_validated],
                                        self.field_values[triangle_without_index[1]],
                                        self.field_values[last_validated],
                                        self.cost(self.density[triangleindex]))
                    if(triangle_without_index[0] in visited_dict): #Update visited dict
                        if((triangleindex in visited_dict[triangle_without_index[0]])
                            and (visited_dict[triangle_without_index[0]][triangleindex]
                                 == self.field_values[triangle_without_index[0]])):
                            self.field_values[triangle_without_index[0]] = min([
                                visited_dict[triangle_without_index[0]][tridex]
                                for tridex in visited_dict[triangle_without_index[0]]
                                if tridex != triangleindex]+[potential_value])
                            visited = self.add_in_order_by_field_value(visited,
                                                                  triangle_without_index[0])
                        visited_dict[triangle_without_index[0]][triangleindex] = potential_value
                    else:
                        visited_dict[triangle_without_index[0]] = {}
                        visited_dict[triangle_without_index[0]][triangleindex] = potential_value

                    if potential_value < self.field_values[triangle_without_index[0]]:
                        self.field_values[triangle_without_index[0]] = potential_value
                        visited = self.add_in_order_by_field_value(visited,
                                                              triangle_without_index[0])

                elif(state_map[triangle_without_index[1]] == 2
                     and state_map[triangle_without_index[0]] == 2):

                    potential_value = (self.field_values[last_validated]
                                       + self.cost(self.density[triangleindex])
                                            *EikoSolver.compute_height_length(self.mesh.vertices[last_validated],
                                                                            self.mesh.vertices[triangle_without_index[0]],
                                                                            self.mesh.vertices[triangle_without_index[1]]))


                    for i in range(2):
                        if( triangle_without_index[i] not in visited_dict):
                            visited_dict[triangle_without_index[i]] = {}

                        visited_dict[triangle_without_index[i]][triangleindex] = potential_value
                        if potential_value < self.field_values[triangle_without_index[i]]:
                            self.field_values[triangle_without_index[i]] = potential_value
                            visited = self.add_in_order_by_field_value(visited,
                                                                  triangle_without_index[i])

    def compute_field_constrained_dep_2(self) -> None:
        """Compute the approximated solution of the eikonal equation and stores the approximation in "fieldValues".

        The method used is the FMT algorithm with constraints on the gradient and
        vertices at distance one or two are considered neighbour.
        i.e. constrained = True and NarrowBandDepth = 2 see "opt"
        """
        state_map = [0 for i in range(len(self.mesh.vertices))]

        visited = []
        validated = []
        self.field_values.set_infinity()

        #Dirichlet free exit conditions
        for index in self.mesh.exit_vertices:
            self.field_values[index] = 0
            validated.append(index)
            state_map[index] = 1

        #Construction of the narrow band
        for index in validated:
            for [triangleindex, triangle_without_index] in self.mesh.triangles_per_vertex[index]:

                if(state_map[triangle_without_index[1]] == 0):
                    state_map[triangle_without_index[1]] = 2

                if(state_map[triangle_without_index[0]] == 0):
                    state_map[triangle_without_index[0]] = 2

                if(state_map[triangle_without_index[0]] == 1
                   and state_map[triangle_without_index[1]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_constrained(
                                        self.mesh.vertices[triangle_without_index[1]],
                                        self.mesh.vertices[triangle_without_index[0]],
                                        self.mesh.vertices[index],
                                        self.field_values[triangle_without_index[0]],
                                        self.field_values[index],
                                        self.cost(self.density[triangleindex]))
                    if potential_value < self.field_values[triangle_without_index[1]]:
                        self.field_values[triangle_without_index[1]] = potential_value
                        visited = self.add_in_order_by_field_value(visited,triangle_without_index[1])


                elif(state_map[triangle_without_index[1]] == 1
                     and state_map[triangle_without_index[0]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_constrained(
                                                self.mesh.vertices[triangle_without_index[0]],
                                                self.mesh.vertices[triangle_without_index[1]],
                                                self.mesh.vertices[index],
                                                self.field_values[triangle_without_index[1]],
                                                self.field_values[index],
                                                self.cost(self.density[triangleindex]))
                    if potential_value < self.field_values[triangle_without_index[0]]:
                        self.field_values[triangle_without_index[0]] = potential_value
                        visited = self.add_in_order_by_field_value(visited,triangle_without_index[0])

                elif(state_map[triangle_without_index[1]] == 2
                     and state_map[triangle_without_index[0]] == 2):
                    for other_index in triangle_without_index:
                        for edge_index in self.mesh.triangles_with_edges[triangleindex]:
                            if( self.mesh.edges[edge_index][0] in [index,other_index]
                                and self.mesh.edges[edge_index][1] in [other_index,index]):
                                potential_value = (self.field_values[index]
                                                   + self.cost(self.density[triangleindex])
                                                        *self.mesh.edge_length[edge_index])
                                break
                        if potential_value < self.field_values[other_index]:
                            self.field_values[other_index] = potential_value
                            visited = self.add_in_order_by_field_value(visited,other_index)


        #Main loop
        last_validated = validated[0]
        num_vertices = len(state_map)

        while(len(validated) < num_vertices):

            #Validation of a vertex
            if(len(visited) == 0):
                if(self.opt["debugging"]):
                    print("Completion : ", len(validated), "/", num_vertices)
                msg = "The fast marching ended too fast... Is the Mesh connected ?"
                raise ValueError(msg)
            if(self.field_values[visited[0]] == float("inf")):
                if(self.opt["debugging"]):
                    print("Infinity value detected at vertex (",
                          self.mesh.vertices[visited[0]][0],",",
                          self.mesh.vertices[visited[0]][1], ")" )
                msg = "Gradient computation impossible."
                raise ValueError(msg)

            validated.append(visited[0])
            last_validated = visited[0]
            visited = visited[1:]
            state_map[last_validated] = 1

            #recomputation of the narrow band
            for [triangleindex, triangle_without_index] in self.mesh.triangles_per_vertex[last_validated]:

                if(state_map[triangle_without_index[1]] == 0):
                    state_map[triangle_without_index[1]] = 2

                if(state_map[triangle_without_index[0]] == 0):
                    state_map[triangle_without_index[0]] = 2

                if(state_map[triangle_without_index[0]] == 1
                   and state_map[triangle_without_index[1]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_constrained(
                                                    self.mesh.vertices[triangle_without_index[1]],
                                                    self.mesh.vertices[triangle_without_index[0]],
                                                    self.mesh.vertices[last_validated],
                                                    self.field_values[triangle_without_index[0]],
                                                    self.field_values[last_validated],
                                                    self.cost(self.density[triangleindex]))
                    if potential_value < self.field_values[triangle_without_index[1]]:
                        self.field_values[triangle_without_index[1]] = potential_value
                        visited = self.add_in_order_by_field_value(visited, triangle_without_index[1])


                elif(state_map[triangle_without_index[1]] == 1
                     and state_map[triangle_without_index[0]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_constrained(
                                                    self.mesh.vertices[triangle_without_index[0]],
                                                    self.mesh.vertices[triangle_without_index[1]],
                                                    self.mesh.vertices[last_validated],
                                                    self.field_values[triangle_without_index[1]],
                                                    self.field_values[last_validated],
                                                    self.cost(self.density[triangleindex]))
                    if potential_value < self.field_values[triangle_without_index[0]]:
                        self.field_values[triangle_without_index[0]] = potential_value
                        visited = self.add_in_order_by_field_value(visited, triangle_without_index[0])

                elif(state_map[triangle_without_index[1]] == 2
                     and state_map[triangle_without_index[0]] == 2):
                    for other_index in triangle_without_index:
                        for edge_index in self.mesh.triangles_with_edges[triangleindex]:
                            if (self.mesh.edges[edge_index][0] in [last_validated,other_index]
                                and self.mesh.edges[edge_index][1] in [other_index,last_validated]):
                                potential_value = (self.field_values[last_validated]
                                                   + self.cost(self.density[triangleindex])
                                                        *self.mesh.edge_length[edge_index])
                                break
                        if potential_value < self.field_values[other_index]:
                            self.field_values[other_index] = potential_value
                            visited = self.add_in_order_by_field_value(visited,other_index)

    def compute_field_unconstrained_dep_1(self) -> None:
        """Compute the approximated solution of the eikonal equation and stores the approximation in "fieldValues".

        The method used is the FMT algorithm with constraints on the gradient and
        vertices at distance one or two are considered neighbour.
        i.e. constrained = False and NarrowBandDepth = 1 see "opt"
        """
        state_map = [0 for i in range(len(self.mesh.vertices))]

        visited = []
        validated = []
        self.field_values.set_infinity()

        #Dirichlet free exit condition
        for index in self.mesh.exit_vertices:
            self.field_values[index] = 0
            validated.append(index)
            state_map[index] = 1

        #Computation of the narrow band
        for index in validated:
            for [triangleindex, triangle_without_index] in self.mesh.triangles_per_vertex[index]:

                if(state_map[triangle_without_index[1]] == 0
                   and state_map[triangle_without_index[0]] == 1):
                    state_map[triangle_without_index[1]] = 2

                if(state_map[triangle_without_index[0]] == 0
                   and state_map[triangle_without_index[1]] == 1):
                    state_map[triangle_without_index[0]] = 2

                if(state_map[triangle_without_index[0]] == 1
                   and state_map[triangle_without_index[1]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_unconstrained(
                                                    self.mesh.vertices[triangle_without_index[1]],
                                                    self.mesh.vertices[triangle_without_index[0]],
                                                    self.mesh.vertices[index],
                                                    self.field_values[triangle_without_index[0]],
                                                    self.field_values[index],
                                                    self.cost(self.density[triangleindex]))
                    if potential_value < self.field_values[triangle_without_index[1]]:
                        self.field_values[triangle_without_index[1]] = potential_value
                        visited = self.add_in_order_by_field_value(visited,triangle_without_index[1])

                elif(state_map[triangle_without_index[1]] == 1
                     and state_map[triangle_without_index[0]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_unconstrained(
                                                    self.mesh.vertices[triangle_without_index[0]],
                                                    self.mesh.vertices[triangle_without_index[1]],
                                                    self.mesh.vertices[index],
                                                    self.field_values[triangle_without_index[1]],
                                                    self.field_values[index],
                                                    self.cost(self.density[triangleindex]))
                    if potential_value < self.field_values[triangle_without_index[0]]:
                        self.field_values[triangle_without_index[0]] = potential_value
                        visited = self.add_in_order_by_field_value(visited,triangle_without_index[0])

        #Main loop
        last_validated = validated[0]
        num_vertices = len(state_map)

        while(len(validated) < num_vertices):

            #Validation of a vertex
            if(len(visited) == 0):
                if(self.opt["debugging"]):
                    print("Completion : ", len(validated), "/", num_vertices)
                msg = "The fast marching ended too fast... Is the Mesh connected ?"
                raise ValueError(msg)
            if(self.field_values[visited[0]] == float("inf")):
                if(self.opt["debugging"]):
                    print("Infinity value detected at vertex (",
                          self.mesh.vertices[visited[0]][0],",",
                          self.mesh.vertices[visited[0]][1], ")" )
                msg = "Gradient computation impossible."
                raise ValueError(msg)

            validated.append(visited[0])
            last_validated = visited[0]
            visited = visited[1:]
            state_map[last_validated] = 1

            #Recomputation of the narrow band
            for [triangleindex, triangle_without_index] in self.mesh.triangles_per_vertex[last_validated]:

                if(state_map[triangle_without_index[1]] == 0
                   and state_map[triangle_without_index[0]] == 1):
                    state_map[triangle_without_index[1]] = 2

                if(state_map[triangle_without_index[0]] == 0
                   and state_map[triangle_without_index[1]] == 1):
                    state_map[triangle_without_index[0]] = 2

                if(state_map[triangle_without_index[0]] == 1
                   and state_map[triangle_without_index[1]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_unconstrained(
                                                    self.mesh.vertices[triangle_without_index[1]],
                                                    self.mesh.vertices[triangle_without_index[0]],
                                                    self.mesh.vertices[last_validated],
                                                    self.field_values[triangle_without_index[0]],
                                                    self.field_values[last_validated],
                                                    self.cost(self.density[triangleindex]))
                    if potential_value < self.field_values[triangle_without_index[1]]:
                        self.field_values[triangle_without_index[1]] = potential_value
                        visited = self.add_in_order_by_field_value(visited, triangle_without_index[1])



                elif(state_map[triangle_without_index[1]] == 1
                     and state_map[triangle_without_index[0]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_unconstrained(
                                                    self.mesh.vertices[triangle_without_index[0]],
                                                    self.mesh.vertices[triangle_without_index[1]],
                                                    self.mesh.vertices[last_validated],
                                                    self.field_values[triangle_without_index[1]],
                                                    self.field_values[last_validated],
                                                    self.cost(self.density[triangleindex]))
                    if potential_value < self.field_values[triangle_without_index[0]]:
                        self.field_values[triangle_without_index[0]] = potential_value
                        visited = self.add_in_order_by_field_value(visited, triangle_without_index[0])


    def compute_field_constrained_dep_1(self) -> None:
        """Compute the approximated solution of the eikonal equation and stores the approximation in "fieldValues".

        The method used is the FMT algorithm with constraints on the gradient and
        vertices at distance one or two are considered neighbour.
        i.e. constrained = True and NarrowBandDepth = 1 see "opt"
        """
        state_map = [0 for i in range(len(self.mesh.vertices))]

        visited = []
        validated = []
        self.field_values.set_infinity()

        #Dirichlet free exit condition
        for index in self.mesh.exit_vertices:
            self.field_values[index] = 0
            validated.append(index)
            state_map[index] = 1

        #Computation of the narrow band
        for index in validated:
            for [triangleindex, triangle_without_index] in self.mesh.triangles_per_vertex[index]:

                if(state_map[triangle_without_index[1]] == 0):
                    state_map[triangle_without_index[1]] = 2

                if(state_map[triangle_without_index[0]] == 0):
                    state_map[triangle_without_index[0]] = 2

                if(state_map[triangle_without_index[0]] == 1
                   and state_map[triangle_without_index[1]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_constrained(
                                                    self.mesh.vertices[triangle_without_index[1]],
                                                    self.mesh.vertices[triangle_without_index[0]],
                                                    self.mesh.vertices[index],
                                                    self.field_values[triangle_without_index[0]],
                                                    self.field_values[index],
                                                    self.cost(self.density[triangleindex]))
                    if potential_value < self.field_values[triangle_without_index[1]]:
                        self.field_values[triangle_without_index[1]] = potential_value
                        visited = self.add_in_order_by_field_value(visited,triangle_without_index[1])


                elif(state_map[triangle_without_index[1]] == 1
                     and state_map[triangle_without_index[0]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_constrained(
                                                    self.mesh.vertices[triangle_without_index[0]],
                                                    self.mesh.vertices[triangle_without_index[1]],
                                                    self.mesh.vertices[index],
                                                    self.field_values[triangle_without_index[1]],
                                                    self.field_values[index],
                                                    self.cost(self.density[triangleindex]))
                    if potential_value < self.field_values[triangle_without_index[0]]:
                        self.field_values[triangle_without_index[0]] = potential_value
                        visited = self.add_in_order_by_field_value(visited,triangle_without_index[0])

                elif(state_map[triangle_without_index[1]] == 2
                     and state_map[triangle_without_index[0]] == 2):
                    potential_value = (self.field_values[index]
                                       + self.cost(self.density[triangleindex])
                                            *EikoSolver.compute_height_length(self.mesh.vertices[index],
                                                                            self.mesh.vertices[triangle_without_index[0]],
                                                                            self.mesh.vertices[triangle_without_index[1]]))

                    for i in range(2):
                        if potential_value < self.field_values[triangle_without_index[i]]:
                            self.field_values[triangle_without_index[i]] = potential_value
                            visited = self.add_in_order_by_field_value(visited, triangle_without_index[i])

        #Main loop
        last_validated = validated[0]
        num_vertices = len(state_map)

        while(len(validated) < num_vertices):

            #Validation of a vertex
            if(len(visited) == 0):
                if(self.opt["debugging"]):
                    print("Completion : ", len(validated), "/", num_vertices)
                msg = "The fast marching ended too fast... Is the Mesh connected ?"
                raise ValueError(msg)
            if(self.field_values[visited[0]] == float("inf")):
                if(self.opt["debugging"]):
                    print("Infinity value detected at vertex (",
                          self.mesh.vertices[visited[0]][0],",",
                          self.mesh.vertices[visited[0]][1], ")" )
                msg = "Gradient computation impossible."
                raise ValueError(msg)

            validated.append(visited[0])
            last_validated = visited[0]
            visited = visited[1:]
            state_map[last_validated] = 1

            #Recomputation of the narrow band
            for [triangleindex, triangle_without_index] in self.mesh.triangles_per_vertex[last_validated]:

                if(state_map[triangle_without_index[1]] == 0):
                    state_map[triangle_without_index[1]] = 2

                if(state_map[triangle_without_index[0]] == 0):
                    state_map[triangle_without_index[0]] = 2

                if(state_map[triangle_without_index[0]] == 1
                   and state_map[triangle_without_index[1]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_constrained(
                                                    self.mesh.vertices[triangle_without_index[1]],
                                                    self.mesh.vertices[triangle_without_index[0]],
                                                    self.mesh.vertices[last_validated],
                                                    self.field_values[triangle_without_index[0]],
                                                    self.field_values[last_validated],
                                                    self.cost(self.density[triangleindex]))
                    if potential_value < self.field_values[triangle_without_index[1]]:
                        self.field_values[triangle_without_index[1]] = potential_value
                        visited = self.add_in_order_by_field_value(visited,
                                                              triangle_without_index[1])


                elif(state_map[triangle_without_index[1]] == 1
                     and state_map[triangle_without_index[0]] == 2):
                    potential_value = EikoSolver.compute_height_from_grad_constrained(
                                                    self.mesh.vertices[triangle_without_index[0]],
                                                    self.mesh.vertices[triangle_without_index[1]],
                                                    self.mesh.vertices[last_validated],
                                                    self.field_values[triangle_without_index[1]],
                                                    self.field_values[last_validated],
                                                    self.cost(self.density[triangleindex]))
                    if potential_value < self.field_values[triangle_without_index[0]]:
                        self.field_values[triangle_without_index[0]] = potential_value
                        visited = self.add_in_order_by_field_value(visited,
                                                              triangle_without_index[0])

                elif(state_map[triangle_without_index[1]] == 2
                     and state_map[triangle_without_index[0]] == 2):
                    potential_value = (self.field_values[last_validated]
                                       + self.cost(self.density[triangleindex])
                                            *EikoSolver.compute_height_length(self.mesh.vertices[last_validated],
                                                                            self.mesh.vertices[triangle_without_index[0]],
                                                                            self.mesh.vertices[triangle_without_index[1]]))

                    for i in range(2):
                        if potential_value < self.field_values[triangle_without_index[i]]:
                            self.field_values[triangle_without_index[i]] = potential_value
                            visited = self.add_in_order_by_field_value(visited,
                                                                  triangle_without_index[i])


    def compute_field_by_edges(self) -> None:
        """Compute a numerical approximation of the solution to the eikonal equation using a FME algorithm."""
        state_map = [0 for i in range(len(self.mesh.vertices))]

        visited = []
        validated = []
        self.field_values.set_infinity()

        #Dirichlet free exit condition
        for index in self.mesh.exit_vertices:
            self.field_values[index] = 0
            validated.append(index)
            state_map[index] = 1

        #Computation of the narrow band
        for index in validated:
            for [triangleindex, triangle_without_index] in self.mesh.triangles_per_vertex[index]:
                for other_index in triangle_without_index:
                    if(state_map[other_index] == 0):
                        state_map[other_index] = 2
                    for edge_index in self.mesh.triangles_with_edges[triangleindex]:
                        if (self.mesh.edges[edge_index][0] in [index,other_index]
                            and self.mesh.edges[edge_index][1] in [other_index,index]):
                            potential_value = (self.field_values[index]
                                               + self.cost(self.density[triangleindex])
                                               *self.mesh.edge_length[edge_index])
                            break
                    if potential_value < self.field_values[other_index]:
                        self.field_values[other_index] = potential_value
                        visited = self.add_in_order_by_field_value(visited,other_index)

        #Main loop
        last_validated = validated[0]
        num_vertices = len(state_map)

        while(len(validated) < num_vertices):

            #Validation of a vertex
            if(len(visited) == 0):
                if(self.opt["debugging"]):
                    print("Completion : ", len(validated), "/", num_vertices)
                msg = "The fast marching ended too fast... Is the Mesh connected ?"
                raise ValueError(msg)
            if(self.field_values[visited[0]] == float("inf")):
                if(self.opt["debugging"]):
                    print("Infinity value detected at vertex (",
                          self.mesh.vertices[visited[0]][0],",",
                          self.mesh.vertices[visited[0]][1], ")" )
                msg = "Gradient computation impossible."
                raise ValueError(msg)

            validated.append(visited[0])
            last_validated = visited[0]
            visited = visited[1:]
            state_map[last_validated] = 1

            #Recomputation of the narrow band
            for [triangleindex, triangle_without_index] in self.mesh.triangles_per_vertex[last_validated]:

                for other_index in triangle_without_index:
                    if(state_map[other_index] == 0):
                        state_map[other_index] = 2
                    for edge_index in self.mesh.triangles_with_edges[triangleindex]:
                        if (self.mesh.edges[edge_index][0] in [last_validated,other_index]
                            and self.mesh.edges[edge_index][1] in [other_index,last_validated]):
                            potential_value = (self.field_values[last_validated]
                                               + self.cost(self.density[triangleindex])
                                                    *self.mesh.edge_length[edge_index])
                            break
                    if potential_value < self.field_values[other_index]:
                        self.field_values[other_index] = potential_value
                        visited = self.add_in_order_by_field_value(visited,other_index)

    def compute_field(self) -> None:
        """Compute the approximated solution of the eikonal equation and stores the approximation in "fieldValues".

        The method used depends on the options set in the "opt" dictionary.
        """
        if(self.opt["method"] == "FMT"):
            if(self.opt["constrained"] and self.opt["NarrowBandDepth"] == 2):
                self.compute_field_constrained_dep_2()
            elif( not self.opt["constrained"] and self.opt["NarrowBandDepth"] == 2):
                self.compute_field_unconstrained_dep_2()
            elif( self.opt["constrained"] and self.opt["NarrowBandDepth"] == 1):
                self.compute_field_constrained_dep_1()
            elif( not self.opt["constrained"] and self.opt["NarrowBandDepth"] == 1):
                self.compute_field_unconstrained_dep_1()
        elif(self.opt["method"] == "FME"):
            self.compute_field_by_edges()

    def show_narrow_band_after_step(self, n:int) -> None:
        """Compute the first `n` vertices of the approximation, then, display the narrow band.

        Args:
            n (int): number of vertices to compute.

        """
        not_visited = list(range(len(self.mesh.vertices)))

        visited = []
        validated = []
        self.field_values.set_infinity()

        #Dirichlet free exit condition
        for index in self.mesh.exit_vertices:
            self.field_values[index] = 0
            validated.append(index)
            not_visited.remove(index)

        #Computation of the narrow band
        for index in validated:
            for triangleindex, triangle in enumerate(self.mesh.triangles):
                if index in triangle:
                    count = 0
                    selected = -1
                    for j,i in enumerate(triangle):
                        if(i in validated):
                            count += 1
                            offset_validated_point = j
                        elif i in not_visited or i in visited:
                            selected = i
                            offset = j

                    if(count == 2 and selected != -1):
                        if(selected not in visited):
                            visited.append(selected)
                            not_visited.remove(selected)

                        if(self.opt["constrained"]):
                            potential_value = EikoSolver.compute_height_from_grad_constrained(
                                                self.mesh.vertices[triangle[offset]],
                                                self.mesh.vertices[triangle[(offset + 1)%3]],
                                                self.mesh.vertices[triangle[(offset + 2)%3]],
                                                self.field_values[triangle[(offset + 1)%3]],
                                                self.field_values[triangle[(offset + 2)%3]],
                                                self.cost(self.density[triangleindex]))
                        else:
                            potential_value = EikoSolver.compute_height_from_grad_unconstrained(
                                                self.mesh.vertices[triangle[offset]],
                                                self.mesh.vertices[triangle[(offset + 1)%3]],
                                                self.mesh.vertices[triangle[(offset + 2)%3]],
                                                self.field_values[triangle[(offset + 1)%3]],
                                                self.field_values[triangle[(offset + 2)%3]],
                                                self.cost(self.density[triangleindex]))
                        self.field_values[selected] = min(potential_value,
                                                          self.field_values[selected])



                    elif(self.opt["NarrowBandDepth"] == 2 and count == 1):

                        for i in triangle:
                            if(i in not_visited):
                                visited.append(i)
                                not_visited.remove(i)

                            if(i not in validated):
                                potential_value = (self.field_values[triangle[offset_validated_point]]
                                                   + self.cost(self.density[triangleindex])
                                                        *EikoSolver.compute_height_length(
                                                                self.mesh.vertices[triangle[offset_validated_point]],
                                                                self.mesh.vertices[triangle[(offset_validated_point + 1)%3]],
                                                                self.mesh.vertices[triangle[(offset_validated_point + 2)%3]]))
                                self.field_values[i] = min(potential_value,
                                                           self.field_values[i])

        #Main loop
        last_validated = validated[0]
        for vide in range(n):

            #Recomputation of the narrow band
            for triangleindex,triangle in enumerate(self.mesh.triangles):
                if(last_validated in triangle):
                    count = 0
                    selected = -1
                    for j,i in enumerate(triangle):
                        if(i in validated):
                            count += 1
                            offset_validated_point = j
                        elif i in not_visited or i in visited:
                            selected = i
                            offset = j

                    if(count == 2 and selected != -1):
                        if(selected not in visited):
                            visited.append(selected)
                            not_visited.remove(selected)

                        if(self.opt["constrained"]):
                            potential_value = EikoSolver.compute_height_from_grad_constrained(
                                                            self.mesh.vertices[triangle[offset]],
                                                            self.mesh.vertices[triangle[(offset + 1)%3]],
                                                            self.mesh.vertices[triangle[(offset + 2)%3]],
                                                            self.field_values[triangle[(offset + 1)%3]],
                                                            self.field_values[triangle[(offset + 2)%3]],
                                                            self.cost(self.density[triangleindex]))
                        else:
                            potential_value = EikoSolver.compute_height_from_grad_unconstrained(
                                                            self.mesh.vertices[triangle[offset]],
                                                            self.mesh.vertices[triangle[(offset + 1)%3]],
                                                            self.mesh.vertices[triangle[(offset + 2)%3]],
                                                            self.field_values[triangle[(offset + 1)%3]],
                                                            self.field_values[triangle[(offset + 2)%3]],
                                                            self.cost(self.density[triangleindex]))
                        self.field_values[selected] = min(potential_value,
                                                          self.field_values[selected])



                    elif(self.opt["NarrowBandDepth"] == 2 and count == 1):

                        for i in triangle:
                            if(i in not_visited):
                                visited.append(i)
                                not_visited.remove(i)

                            if(i not in validated):
                                potential_value = (self.field_values[triangle[offset_validated_point]]
                                                   + self.cost(self.density[triangleindex])
                                                        *EikoSolver.compute_height_length(
                                                            self.mesh.vertices[triangle[offset_validated_point]],
                                                            self.mesh.vertices[triangle[(offset_validated_point + 1)%3]],
                                                            self.mesh.vertices[triangle[(offset_validated_point + 2)%3]]))
                                self.field_values[i] = min(potential_value,
                                                           self.field_values[i])

            #Validation:
            min_index = visited[0]
            min_value = self.field_values[min_index]
            for index in visited:
                if min_value > self.field_values[index]:
                    min_index = index
                    min_value = self.field_values[index]

            if(min_value == float("inf")):
                if(self.opt["debugging"]):
                    print("Infinity value detected at vertex (",
                          self.mesh.vertices[visited[0]][0],",",
                          self.mesh.vertices[visited[0]][1], ")" )
                msg = "Gradient computation impossible."
                raise ValueError(msg)

            validated.append(min_index)
            last_validated = min_index
            visited.remove(min_index)

        visited.append(min_index)
        validated.remove(min_index)

        for index in validated:
            self.field_values[index] = 0.5
        for index in visited:
            self.field_values[index] = 0
        for index in not_visited:
            self.field_values[index] = 1

        self.field_values.show(grid=True)

    def add_in_order_by_field_value(self,ordered_list:list[int],index:int) -> list[int]:
        """Add the index given as an argument to the given list.

        Returns the new list. The index is inserted such that the `fieldValues`
        corresponding to the vertices of the indices are increasing.

        Example:
            If we have the following lists::

                MyMesh = Mesh()
                MyMesh.vertices = [[0,0],[0,1],[1,1],[1,0]]
                MySolver = EikoSolver(MyMesh)
                MySolver.fieldValues.values = [0.5,4.5,2.0,3.0]
                Mylist = [0,3]

            Then::

                >>>MySolver.add_in_order_by_field_value(Mylist,1)
                [0,3,1]
                >>>MySolver.add_in_order_by_field_value(Mylist,2)
                [0,2,3]

        Note:
            If an index is already present in the list, its position in the list is
            recomputed.

        Args:
            ordered_list (list[int]): the list in which the new index will be inserted.
            index (int): the index to insert.

        Returns:
            list[int]: the new list containing the index.

        """
        try: #Kept as opposed to SIM105 recommandation for performance reasons.
            ordered_list.remove(index)
        except ValueError:
            pass

        if(len(ordered_list) == 0):
            return [index]
        rank = 0
        while(self.field_values[index] > self.field_values[ordered_list[rank]]):
            rank += 1
            if rank >= len(ordered_list):
                break

        return [*ordered_list[:rank],
                index,
                *ordered_list[rank:]]


    @staticmethod
    def compute_height_from_grad_unconstrained(c:list[float],b:list[float],a:list[float],v_b:float,v_a:float,slope:float) -> float:
        r"""Compute and return the value of v_c.

        v_c is computed such that, if we denote by :math:`\\Phi_{ABC}(v_a,v_b,v_c)`
        the unique affine function :math:`F` of :math:`\\mathbb{R}^2` such that

        .. math::
            F(A) = v_a, \\; F(B) = v_b, \\; F(C) = v_c

        then we have

        .. math::
            \\left| \\nabla\\Phi_{ABC}(v_a,v_b,v_c) \\right| = slope.

        Args:
            c (list[float]): a point of the triangle.
            b (list[float]): a point of the triangle.
            a (list[float]): a point of the triangle.
            v_b (float): the value of F at b.
            v_a (float): the value of F at a.
            slope (float): the norm of the gradient.

        Returns:
            float: the value v_c.

        """
        if(v_b == float("inf") or v_a == float("inf") or slope < 0):
            msg = (f"Invalid values for the gradient constraint:"
                   f" v_b = {v_b}, v_a = {v_a}, slope = {slope}")
            raise ValueError(msg)
        ab2 = (a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1])


        if(ab2*(slope**2) - ( (v_a - v_b)**2) <0):
            msg = (f"Warning - invalid values for the gradient constraint:"
                   f" v_b = {v_b}, v_a = {v_a}, slope = {slope}")
            print(msg)
            ac2 = (a[0] - c[0])*(a[0] - c[0]) + (a[1] - c[1])*(a[1] - c[1])
            return v_a + np.sqrt(ac2)*slope


        detcacb = abs((c[0] - a[0])*(c[1] - b[1]) - (c[1] - a[1])*(c[0] - b[0]))
        cbxca = (b[0] - c[0])*(a[0] - c[0]) + (b[1] - c[1])*(a[1] - c[1])

        if(v_a > v_b):
            bc2 = (c[0] - b[0])*(c[0] - b[0]) + (c[1] - b[1])*(c[1] - b[1])
            return (v_b + (v_a-v_b)*(bc2 - cbxca)/ab2
                    + detcacb*np.sqrt(ab2*(slope**2) - ( (v_a - v_b)**2) )/ab2)
        ac2 = (a[0] - c[0])*(a[0] - c[0]) + (a[1] - c[1])*(a[1] - c[1])
        return (v_a + (v_b-v_a)*(ac2 - cbxca)/ab2
                + detcacb*np.sqrt(ab2*(slope**2) - ( (v_b - v_a)**2) )/ab2)

    @staticmethod
    def compute_height_from_grad_constrained(a:list[float],b:list[float],c:list[float],v_b:float,v_c:float,slope:float) -> float:
        r"""Compute and return the value of v_a.

        v_a is computed such that, if we denote by :math:`\\Phi_{ABC}(v_a,v_b,v_c)`
        the unique affine function :math:`F` of :math:`\\mathbb{R}^2` such that

        .. math::
            F(A) = v_a, \\; F(B) = v_b, \\; F(C) = v_c

        then we have

        .. math::
            \\left| \\nabla\\Phi_{ABC}(v_a,v_b,v_c) \\right| = slope.

        If the gradient found is not inside the triangle, v_a is taken as the smallest
        value found by following an edge of the triangle.

        Args:
            a (list[float]): a point of the triangle.
            b (list[float]): a point of the triangle.
            c (list[float]): a point of the triangle.
            v_b (float): the value of F at b.
            v_c (float): the value of F at c.
            slope (float): the norm of the gradient.

        Returns:
            float: the value v_a.

        """
        if(v_b == float("inf") or v_c == float("inf") or slope < 0):
            msg = (f"Invalid values for the gradient constraint:"
                   f" v_b = {v_b}, v_c = {v_c}, slope = {slope}")
            raise ValueError(msg)

        if( v_b > v_c):
            s1 = v_c
            s2 = list(c)
            c = b
            v_c = v_b
            b = s2
            v_b = s1

        ab2 = (a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1])
        abxbc = (b[0] - c[0])*(a[0] - b[0]) + (b[1] - c[1])*(a[1] - b[1])

        if(slope*abxbc > (v_c-v_b)*np.sqrt(ab2)):
            return v_b + np.sqrt(ab2)*slope

        ac2 = (a[0] - c[0])*(a[0] - c[0]) + (a[1] - c[1])*(a[1] - c[1])
        cbxca = (b[0] - c[0])*(a[0] - c[0]) + (b[1] - c[1])*(a[1] - c[1])

        if(slope*cbxca < (v_c-v_b)*np.sqrt(ac2)):
            return v_c + np.sqrt(ac2)*slope

        bc2 = (b[0] - c[0])*(b[0] - c[0]) + (b[1] - c[1])*(b[1] - c[1])

        if(bc2*slope**2 - ( (v_c - v_b)**2) <0):
            msg = (f"Warning - invalid values for the gradient constraint:"
                   f" v_b = {v_b}, v_c = {v_c}, slope = {slope}")
            print(msg)
            return v_c + np.sqrt(ac2)*slope

        detabcb = abs((b[0] - a[0])*(b[1] - c[1]) - (b[1] - a[1])*(b[0] - c[0]))

        return (v_b - (v_c-v_b)*(abxbc)/bc2
                + detabcb*np.sqrt(bc2*slope**2 - ( (v_b - v_c)**2) )/bc2)

    @staticmethod
    def compute_height_length(b:float,c:float,a:float) -> float:
        """Compute and return the length of the height of the triangle ABC passing through B.

        Args:
            b (list[float]): a point of the triangle.
            c (list[float]): a point of the triangle.
            a (list[float]): a point of the triangle.

        Returns:
            float: the length of the height.

        """
        detabac = abs((a[0] - c[0])*(a[1] - b[1]) - (a[1] - c[1])*(a[0] - b[0]))
        ac2 = (a[0] - c[0])*(a[0] - c[0]) + (a[1] - c[1])*(a[1] - c[1])
        return detabac/np.sqrt(ac2)

    @staticmethod
    def dist(x0:float,y0:float,x1:float,y1:float) -> float:
        """Compute and return the euclidian distance between (x0,y0) and (x1,y1)."""
        return(np.sqrt((x1-x0)**2 + (y1-y0)**2))
