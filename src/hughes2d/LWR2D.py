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
from functools import partial

import numpy as np
from numpy.typing import ArrayLike

from hughes2d.Mesh2D import CellValueMap, Mesh

PRECISION = 1e-10
DEFAULT_PRECISION = 1e-5
EMPTY_LIST = []
DEFAULT_OPTIONS = {"convexFlux" : True, "anNum" : "dichotomy", "method" : "midVector"}

class LWRSolver:
    r"""A solver for LWR-type scalar conservation law directed by a vector field in two dimensions.

    We use the following notations:

    .. math::

        \left\{ \begin{matrix}
        \partial_t \rho(t,x) + \mathbf{div}\left[ \rho(t,x) v(\rho(t,x)) \vec{V}(t,x) \right] = 0, \\
        \rho(0,x) = \rho_0(x)
        \end{matrix} \right.

    Args:
        Mesh (Mesh): a mesh object on which the equation will be approximated.
        previous_density (CellValueMap or list[float]): initial density for the solver.
            Must be of the shape of `Mesh.triangles`. Represents :math:`\rho_0(x)` in
                the equation above.
        direction_map (list[list[float]]): a vector field represented by a list of
            vectors with the shape of `Mesh.triangles`. Typically this corresponds to
            the output of `VertexValueMap.computeGradientFlow()`.
            Represents :math:`\vec{V}(t,x)` in the equation above.
        speed_function (function, float -> float): a function respresenting the speed of
            the agents depending on the local density. Represents :math:`v(\cdot)` in
            the equation above.
        dt (float): the time division for the approximation.

            Warning:
                The CFL condition must be satisfied for the simulations to make sense.
                Here the CFL condition is:

                .. math::
                    \Delta t \leq \frac{\underline{|\triangle|}}{3\underline{\textrm{$\triangle$}}Lip_f},

                where :math:`\underline{|\triangle|}` denotes the minimal area of a
                triangle in the mesh :math:`M_\Delta` and
                :math:`\underline{\textrm{$\triangle$}}` denotes the maximal length
                of the edges of the mesh.

        opt (dict): an optional dictionary prescribing the numerical method.

            ========================= ====== ====================== =========================================================================================
            key                       type   possible values        description
            ========================= ====== ====================== =========================================================================================
            "method"                  str    "tmap" or "midvector"  determines of the conflict between non-colinear vectors are resolved
            "anNum"                   str    "dichotomy"            only parameter available at the moment, numerical method to use for the approximations
            "convexFlux"              bool   True or False          optimization of the computations when the flux is convex or concave
            "ApproximationThreshold"  float  > 0                    the precision to use for the numerical approximations
            "debugging"               bool   True or False          determines if the solver will print debugging informations in the console or not
            ========================= ====== ====================== =========================================================================================

    Attributes
    ----------

    Attributes:
        mesh (Mesh): a mesh object on which the equation will be approximated.
        densityt0 (ArrayLike): initial density for the solver. Of the shape of
            `Mesh.triangles`. Represents :math:`\rho_0(x)` in the equation above.
        densityt1 (ArrayLike): density computed by the solver after one time step.
            Of the shape of `Mesh.triangles`. Represents :math:`\rho(\Delta t, x)`.
        directions (list[list[float]]): a vector field represented by a list of vectors
            with the shape of `Mesh.triangles`. Typically this corresponds to the output
            of `VertexValueMap.computeGradientFlow()`. Represents :math:`\vec{V}(t,x)`
            in the equation above.
        fluxFunction (function, float -> float): a function respresenting the flux of
            agents depending on the local density. Represents
            :math:`\rho \mapsto \rho v(\cdot)` in the equation above.
        dt (float): the time division for the approximation.

            Warning:
                The CFL condition must be satisfied for the simulations to make sense.
                Here the CFL condition is:

                .. math::
                    \Delta t \leq \frac{\underline{|\triangle|}}{3\underline{\textrm{$\triangle$}}Lip_f},

                where :math:`\underline{|\triangle|}` denotes the minimal area of a
                triangle in the mesh :math:`M_\Delta` and
                :math:`\underline{\textrm{$\triangle$}}` denotes the maximal length
                of the edges of the mesh.

        opt (dict): an optional dictionary prescribing the numerical method.

    Methods
    -------
    """

    def __init__(self,  mesh:Mesh,
                        dt:float,
                        previous_density: list[float] = EMPTY_LIST,
                        direction_map:list[list[float]] = EMPTY_LIST,
                        speed_function: Callable[[float], float] = (lambda x: 1-x),
                        opt:dict = DEFAULT_OPTIONS) -> None:
        self.mesh:Mesh = mesh
        if(previous_density != EMPTY_LIST):
            self.densityt0:ArrayLike = np.array(previous_density) #type CellValueMap
        else:
            self.densityt0:ArrayLike = np.array([0.0 for t in self.mesh.triangles])
        self.densityt1:ArrayLike = np.array([0.0 for t in self.mesh.triangles])

        self.dt:float = dt

        self.opt:dict = opt

        self.directions:list[list[float]] = direction_map

        if("ApproximationThreshold" not in self.opt ):
            self.opt["ApproximationThreshold"] = DEFAULT_PRECISION

        if("debugging" not in self.opt ):
            self.opt["debugging"] = False

        self.fluxFunction = lambda x: x*speed_function(x)

        if(self.opt["convexFlux"]):
            if(self.fluxFunction == (lambda x: x*(1-x))):
                self.max_flux_point = 0.5
                if self.opt["debugging"]:
                    print("Explicit max point considered to be 0.5")
            else:
                self.max_flux_point = LWRSolver.arg_max(self.fluxFunction,0,1,
                                                     self.opt["ApproximationThreshold"])
                if self.opt["debugging"]:
                    print(f"Maximal flux point found at p = {self.max_flux_point}")

    def check_cfl(self, lip_constant:float = 1) -> bool:
        r"""Check if the CFL condition is satisfied.

        Here the CFL condition is:

        .. math::
            \Delta t \leq \frac{\underline{|\triangle|}}{3\underline{\textrm{$\triangle$}}Lip_f},

        where :math:`\underline{|\triangle|}` denotes the minimal area of a triangle
        in the mesh :math:`M_\Delta` and :math:`\underline{\textrm{$\triangle$}}`
        denotes the maximal length of the edges of the mesh.

        Args:
            lip_constant (float): the Lipschitz constant of the flux.

        Returns:
            bool: True if the CFL condition is verified.

        """
        if 3*self.dt*lip_constant*self.mesh.max_edge_length > self.mesh.min_cell_area:
            print("Warning - the CFL conditions might not be verified everywhere in"
                  " the mesh.")
            if self.opt["debugging"]:
                print(f"CFL values : dt = {self.dt} VS"
                      f"{self.mesh.min_cell_area/(3*lip_constant*self.mesh.max_edge_length)}")
            return False
        return True


    def compute_step_tmap(self) -> None:
        """Compute the approximated solution to the scalar conservation law after one time step with the transmission maps method."""
        if(self.opt["convexFlux"]):
            for triangle_index, triangle_cell in enumerate(self.mesh.triangles_with_edges):
                modif_density = 0
                for edge_number, edge_index in enumerate(triangle_cell):
                    pair_of_triangles = self.mesh.pairs_of_triangles[edge_index]
                    triangle_grad = self.directions[triangle_index]
                    normal = self.mesh.outer_normal_vect_by_triangles[triangle_index][edge_number]

                    triangle_flux = (triangle_grad[0]*normal[0]
                                     + triangle_grad[1]*normal[1])

                    if(edge_index in self.mesh.exit_edges):
                        far_flux = 1
                        far_density = 0
                    elif(edge_index in self.mesh.wall_edges):
                        triangle_flux = 0
                        far_flux = 0
                        far_density = 0
                    else:
                        if(pair_of_triangles[0] == triangle_index):
                            far_triangle_grad = self.directions[pair_of_triangles[1]]
                            far_density = self.densityt0[pair_of_triangles[1]]
                        else:
                            far_triangle_grad = self.directions[pair_of_triangles[0]]
                            far_density = self.densityt0[pair_of_triangles[0]]
                        far_flux = (far_triangle_grad[0]*normal[0]
                                    + far_triangle_grad[1]*normal[1])

                    if (triangle_flux > 1 + self.opt["ApproximationThreshold"]
                        or far_flux > 1 + self.opt["ApproximationThreshold"]):
                        print("Warning : scalar product out of bonds")
                        if self.opt["debugging"]:
                            print("| far_triangle_grad | = "
                                  f"{(far_triangle_grad[0]*far_triangle_grad[0]+ far_triangle_grad[1]*far_triangle_grad[1])}")
                            print(f"| normal | = {normal[0]*normal[0]+normal[1]*normal[1]}")


                    """
                    We search for k such that
                    god(triangle_flux function, densityt0[triangle_index], k )
                    = god(far_flux function, k, fardensity).
                    In the specific case of a convex flux function we can solve that by
                    treating different cases :
                    """
                    if((triangle_flux <= 0 and far_flux >= 0)
                        or (triangle_flux >= 0 and far_flux <= 0)):
                        total_flux = 0

                    elif(triangle_flux < 0 and far_flux < 0):
                        outing_flux = (triangle_flux
                                       *self.fluxFunction(max(self.densityt0[triangle_index],
                                                              self.max_flux_point)))
                        entering_flux = far_flux*self.fluxFunction(min(far_density,
                                                                       self.max_flux_point))
                        total_flux = max(outing_flux,entering_flux)

                    elif(triangle_flux > 0 and far_flux > 0):
                        entering_flux = (far_flux*self.fluxFunction(max(far_density,
                                                                        self.max_flux_point)))
                        outing_flux = (triangle_flux
                                       *self.fluxFunction(min(self.densityt0[triangle_index],
                                                              self.max_flux_point)))
                        total_flux = min(outing_flux,entering_flux)
                    modif_density -= ((self.dt/self.mesh.cell_areas[triangle_index])
                                      *self.mesh.edge_length[edge_index]*total_flux)
                self.densityt1[triangle_index] = (self.densityt0[triangle_index]
                                                  +modif_density)

        else:
            for triangle_index, triangle_cell in enumerate(self.mesh.triangles_with_edges):
                modif_density = 0
                for edge_number, edge_index in enumerate(triangle_cell):
                    pair_of_triangles = self.mesh.pairs_of_triangles[edge_index]
                    triangle_grad = self.directions[triangle_index]
                    normal = self.mesh.outer_normal_vect_by_triangles[triangle_index][edge_number]

                    triangle_flux = triangle_grad[0]*normal[0] + triangle_grad[1]*normal[1]

                    if(edge_index in self.mesh.exit_edges):
                        far_flux = 1
                        far_density = 0
                    elif(edge_index in self.mesh.wall_edges):
                        triangle_flux = 0
                        far_flux = 0
                        far_density = 0
                    else:
                        if(pair_of_triangles[0] == triangle_index):
                            far_triangle_grad = self.directions[pair_of_triangles[1]]
                            far_density = self.densityt0[pair_of_triangles[1]]
                        else:
                            far_triangle_grad = self.directions[pair_of_triangles[0]]
                            far_density = self.densityt0[pair_of_triangles[0]]
                        far_flux = (far_triangle_grad[0]*normal[0] + far_triangle_grad[1]*normal[1])

                    if (triangle_flux > 1 + self.opt["ApproximationThreshold"]
                        or far_flux > 1 + self.opt["ApproximationThreshold"]):
                        print("Warning : scalar product out of bonds")
                        if self.opt["debugging"]:
                            print(f"| far_triangle_grad | = {far_triangle_grad[0]*far_triangle_grad[0]+far_triangle_grad[1]*far_triangle_grad[1]}")
                            print(f"| normal | = {normal[0]*normal[0]+normal[1]*normal[1]}")
                    """
                    If the flux function is not convex, we use a dichotomy to solve
                    the problem.
                    """
                    if((triangle_flux <= 0 and far_flux >= 0)
                        or (triangle_flux >= 0 and far_flux <= 0)):
                        total_flux = 0
                    else:
                        def linear_func(x:float, a:float) -> float:
                            return a*self.fluxFunction(x)
                        outer_flux = partial(linear_func, a=triangle_flux)
                        inner_flux = partial(linear_func, a=far_flux)

                        def transmission_map(k:float,
                                             density_1:float,
                                             flux_function_1:Callable[[float],float],
                                             density_2:float,
                                             flux_function_2:Callable[[float],float]) -> float:
                            return (LWRSolver.god(flux_function_1, density_1, k,
                                                 precision=self.opt["ApproximationThreshold"])
                                 - LWRSolver.god(flux_function_2,k, density_2,
                                                 precision=self.opt["ApproximationThreshold"]))

                        parametrized_flux = partial(transmission_map,
                                                    density_1 = self.densityt0[triangle_index],
                                                    flux_function_1 = outer_flux,
                                                    density_2 = far_density,
                                                    flux_function_2 = inner_flux)
                        k = LWRSolver.appro_zero_dichotomy(parametrized_flux,0,1,
                                                          self.opt["ApproximationThreshold"],
                                                          hints=[self.densityt0[triangle_index],far_density])
                        total_flux = LWRSolver.god(outer_flux,
                                                   self.densityt0[triangle_index],
                                                   k,
                                                   self.opt["ApproximationThreshold"])
                    modif_density -= ((self.dt/self.mesh.cell_areas[triangle_index])
                                        *self.mesh.edge_length[edge_index]* total_flux)

                self.densityt1[triangle_index] =min(1,max(self.densityt0[triangle_index]
                                                           + modif_density,0))

                if(np.abs(min(1,max(self.densityt0[triangle_index] + modif_density,0))
                          - (self.densityt0[triangle_index] + modif_density))
                   > PRECISION):
                    print("Warning : the new density computed is not between 0 and 1.")
                    if(self.opt["debugging"]):
                        print(f"Computed density: "
                              f"{self.densityt0[triangle_index] + modif_density}")

    def compute_step_mid_vector(self) -> None:
        """Compute the approximated solution to the scalar conservation law after one time step with the mid vector method."""
        if(self.opt["convexFlux"]):
            for triangle_index, triangle_cell in enumerate(self.mesh.triangles_with_edges):
                modif_density = 0

                for edge_number, edge_index in enumerate(triangle_cell):

                    pair_of_triangles = self.mesh.pairs_of_triangles[edge_index]
                    triangle_grad = self.directions[triangle_index]
                    normal = self.mesh.outer_normal_vect_by_triangles[triangle_index][edge_number]
                    total_flux = 1

                    if(edge_index in self.mesh.exit_edges):
                        vector_flux = triangle_grad
                        far_density = 0
                    elif(edge_index in self.mesh.wall_edges):
                        vector_flux = triangle_grad
                        total_flux = 0
                    else:

                        if(pair_of_triangles[0] == triangle_index):
                            far_triangle_grad = self.directions[pair_of_triangles[1]]
                            far_density = self.densityt0[pair_of_triangles[1]]
                        else:
                            far_triangle_grad = self.directions[pair_of_triangles[0]]
                            far_density = self.densityt0[pair_of_triangles[0]]

                        vector_flux = [(triangle_grad[0]*self.densityt0[triangle_index]
                                        + far_triangle_grad[0]*far_density),
                                       (triangle_grad[1]*self.densityt0[triangle_index]
                                        + far_triangle_grad[1]*far_density) ]
                    norm_vector_flux = np.sqrt(vector_flux[0]*vector_flux[0]
                                               + vector_flux[1]*vector_flux[1])

                    if(norm_vector_flux < self.opt["ApproximationThreshold"]
                       or total_flux == 0):
                        total_flux = 0
                    else:
                        scalar_product = (vector_flux[0]*normal[0]
                                          + vector_flux[1]*normal[1])
                        if(scalar_product/norm_vector_flux > 1 + PRECISION):
                            print("Warning : scalar product greater than the norm,"
                                  " perhaps the normal vector is not normalized")
                        if(scalar_product > 0):
                            if(self.densityt0[triangle_index] <= far_density):
                                total_flux = (scalar_product
                                              *min(self.fluxFunction(self.densityt0[triangle_index]),
                                                   self.fluxFunction(far_density))
                                              /norm_vector_flux)
                            elif(far_density < self.max_flux_point
                                 and self.max_flux_point < self.densityt0[triangle_index]):
                                total_flux = (scalar_product
                                              *self.fluxFunction(self.max_flux_point)
                                              /norm_vector_flux)
                            else:
                                total_flux = (scalar_product
                                              *max(self.fluxFunction(self.densityt0[triangle_index]),
                                                   self.fluxFunction(far_density))
                                              /norm_vector_flux)
                        elif(self.densityt0[triangle_index] >= far_density):
                            total_flux = (scalar_product
                                          *min(self.fluxFunction(self.densityt0[triangle_index]),
                                               self.fluxFunction(far_density))
                                          /norm_vector_flux)
                        elif(far_density > self.max_flux_point
                             and self.max_flux_point > self.densityt0[triangle_index]):
                            total_flux = (scalar_product
                                          *self.fluxFunction(self.max_flux_point)
                                          /norm_vector_flux)
                        else:
                            total_flux = (scalar_product
                                          *max(self.fluxFunction(self.densityt0[triangle_index]),
                                               self.fluxFunction(far_density))
                                          /norm_vector_flux)

                    modif_density -= ((self.dt/self.mesh.cell_areas[triangle_index])
                                        *self.mesh.edge_length[edge_index]*total_flux)

                self.densityt1[triangle_index] = min(1,max(self.densityt0[triangle_index]
                                                           + modif_density,0))

                if(np.abs(min(1,max(self.densityt0[triangle_index] + modif_density,0))
                          - (self.densityt0[triangle_index] + modif_density))
                   > PRECISION):
                    print("Warning : the new density computed is not between 0 and 1.")
                    if(self.opt["debugging"]):
                        print(f"Computed density: "
                              f"{self.densityt0[triangle_index] + modif_density}")

        else:
            for triangle_index, triangle_cell in enumerate(self.mesh.triangles_with_edges):
                modif_density = 0
                for edge_number, edge_index in enumerate(triangle_cell):
                    pair_of_triangles = self.mesh.pairs_of_triangles[edge_index]
                    triangle_grad = self.directions[triangle_index]
                    normal = self.mesh.outer_normal_vect_by_triangles[triangle_index][edge_number]
                    total_flux = 1

                    if(edge_index in self.mesh.exit_edges):
                        vector_flux = triangle_grad
                        far_density = 0
                    elif(edge_index in self.mesh.wall_edges):
                        vector_flux = triangle_grad
                        total_flux = 0
                    else:

                        if(pair_of_triangles[0] == triangle_index):
                            far_triangle_grad = self.directions[pair_of_triangles[1]]
                            far_density = self.densityt0[pair_of_triangles[1]]
                        else:
                            far_triangle_grad = self.directions[pair_of_triangles[0]]
                            far_density = self.densityt0[pair_of_triangles[0]]

                        vector_flux = [(triangle_grad[0]*self.densityt0[triangle_index]
                                        + far_triangle_grad[0]*far_density),
                                       (triangle_grad[1]*self.densityt0[triangle_index]
                                        + far_triangle_grad[1]*far_density) ]
                    norm_vector_flux = np.sqrt(vector_flux[0]*vector_flux[0]
                                               + vector_flux[1]*vector_flux[1])

                    if(norm_vector_flux < self.opt["ApproximationThreshold"]
                       or total_flux == 0):
                        total_flux = 0
                    else:
                        if(((vector_flux[0]*normal[0] + vector_flux[1]*normal[1])
                            /norm_vector_flux) > 1 + PRECISION):
                            print("Warning : scalar product greater than the norm,"
                                  " perhaps the normal vector is not normalized")
                        flux_func = partial(lambda x,scalar,norm : scalar/norm*self.fluxFunction(x),
                                            scalar = vector_flux[0]*normal[0] + vector_flux[1]*normal[1],
                                            norm = norm_vector_flux)
                        total_flux = LWRSolver.god(flux_func,
                                                   self.densityt0[triangle_index],
                                                   far_density,
                                                   precision=self.opt["ApproximationThreshold"])

                    modif_density -= ((self.dt/self.mesh.cell_areas[triangle_index])
                                      * self.mesh.edge_length[edge_index]* total_flux)
                self.densityt1[triangle_index] =min(1,max(self.densityt0[triangle_index]
                                                           + modif_density,0))

                if(np.abs(min(1,max(self.densityt0[triangle_index] + modif_density,0))
                          - (self.densityt0[triangle_index] + modif_density))
                   > PRECISION):
                    print("Warning : the new density computed is not between 0 and 1.")
                    if(self.opt["debugging"]):
                        print(f"Computed density: "
                              f"{self.densityt0[triangle_index] + modif_density}")

    def compute_next_step(self) -> None:
        """Compute the approximated solution to the scalar conservation law after one time step and stores it in `densityt1`."""
        if(self.opt["method"] == "tmap"):
            self.compute_step_tmap()
        elif(self.opt["method"] == "midVector"):
            self.compute_step_mid_vector()

    def update(self, new_direction_field:list[list[float]]) -> None:
        r"""Update the direction vector field with the one passed as a parameter and set the actual approximated solution `densityt1` as the initial datum `densityt0`.

        Args:
            new_direction_field (list[list[float]]): a vector field represented by a
                list of vectors with the shape of `Mesh.triangles`. Typically this
                corresponds to the output of `VertexValueMap.computeGradientFlow()`.
                Represents :math:`\vec{V}(t,x)` in the scalar conservation law
                equation.
        """
        self.densityt0 = np.copy(self.densityt1)
        self.directions = new_direction_field

    def show_density(self, t:int=1, preference:str ="plotly") -> None:
        """Display the density map. If t = 0, it shows the initial datum else if t = 1 (default), it shows the density computed after one time step.

        Args:
            t (int, optional): the time at which the density should be dispayed
                (0 or 1).
            preference (str, optional): set to "plotly" or "matplotlib" to chose the
                preferred plotting package. If only one package is installed the
                preference is ignored.

        Raises:
            ValueError: if t is not 0 or 1.

        """
        cell_map = CellValueMap(self.mesh)
        if t == 0:
            cell_map.values = self.densityt0
        elif t == 1:
            cell_map.values = self.densityt1
        else:
            msg = f"Invalid value for t: {t}"
            raise ValueError(msg)

        cell_map.show(preference = preference)

    @staticmethod
    def arg_max(f:Callable[[float],float],a:float,b:float, precision:float = DEFAULT_PRECISION) -> float:
        """Numerically approximate the maximal point of the function `f` between `a` and `b` and returns the maximal argument.

        Args:
            f (function, float -> float): the function `f` for which the armax will be
                computed.
            a (float): one of the bounds for the search domain of the argmax.
            b (float): one of the bounds for the search domain of the argmax.
            precision (float): the error margin for the approximation.

        Returns:
            float: the argmax computed.

        """
        num_slices = int(1+1/precision)
        max_value = -float("inf")
        step = abs(b-a)/num_slices
        xmax = min(a,b)

        for i in range(num_slices):
            test_value = f(min(a,b) + i*step)
            if(test_value > max_value):
                max_value = test_value
                xmax = min(a,b) + i*step
        return(xmax)

    @staticmethod
    def approxi_max(f:Callable[[float],float],a:float,b:float, precision:float = DEFAULT_PRECISION) -> float:
        """Numerically approximates the maximum of the function `f` between `a` and `b` and returns the maximal value.

        Args:
            f (function, float -> float): the function `f` for which the max will be
                computed.
            a (float): one of the bounds for the search domain of the max.
            b (float): one of the bounds for the search domain of the max.
            precision (float): the error margin for the approximation.

        Returns:
            float: the maximal value computed.

        """
        max_value = -float("inf")
        nb_steps:int = int(np.ceil(np.abs(b-a)/precision))
        if nb_steps == 0:
            print("Warning: precision is less than the difference between the bounds of"
                  " the domain")
            return (a+b)/2

        step = np.abs(b-a)/nb_steps

        for i in range(nb_steps):
            max_value = max(max_value, f(min(a,b) + i*step))
        return(max_value)

    @staticmethod
    def approxi_min(f:Callable[[float],float],a:float,b:float, precision:float = DEFAULT_PRECISION) -> float:
        """Numerically approximate the minimum of the function `f` between `a` and `b` and returns the minimal value.

        Args:
            f (function, float -> float): the function `f` for which the min will be
                computed.
            a (float): one of the bounds for the search domain of the min.
            b (float): one of the bounds for the search domain of the min.
            precision (float): the error margin for the approximation.

        Returns:
            float: the minimal value computed.

        """
        min_value = float("inf")
        nb_steps:int = int(np.ceil(np.abs(b-a)/precision))
        if nb_steps == 0:
            print("Warning: precision is less than the difference between the bounds of"
                  " the domain")
            return (a+b)/2

        step = np.abs(b-a)/nb_steps

        for i in range(nb_steps):
            min_value = min(min_value,  f(min(a,b) + i*step))
        return(min_value)

    @staticmethod
    def god(f:Callable[[float],float],a:float,b:float, precision:float = DEFAULT_PRECISION) -> float:
        r"""Numerically approximate the Godunov flux of the function `f` between `a` and `b` defined by the formula below.

        .. math::

            \mathbf{God}_{f}(a,b) = \left\{ \begin{matrix} \mathbf{min}_{c \in [a,b]} f(c) \textrm{ if } a < b \\
            \mathbf{max}_{c \in [b,a]} f(c) \textrm{ if } a > b. \end{matrix} \right.

        Args:
            f (function, float -> float): the function `f` for which the Godunov flux
                will be computed.
            a (float): one of the parameters for the Godunov flux.
            b (float): one of the parameters for the Godunov flux.
            precision (float): the error margin for the approximation.

        Returns:
            float: the Godunov flux computed.

        """
        if np.abs(a - b) < precision:
            return f((a+b)/2)
        if(a < b):
            return LWRSolver.approxi_min(f,a,b, precision)
        return LWRSolver.approxi_max(f,b,a, precision)

    @staticmethod
    def appro_zero_dichotomy(f:Callable[[float],float],a:float,b:float, precision:float = DEFAULT_PRECISION, hints:list[float] = EMPTY_LIST) -> float:
        """Numerically approximate the a root of the function `f` between `a` and `b` using a dichotomy method.

        Args:
            f (function, float -> float): the function `f` for which a root will be
                computed.
            a (float): one of the bounds for the search domain of the root.
            b (float): one of the bounds for the search domain of the root.
            precision (float,optional): the error margin for the approximation.
            hints (list[float],optional): possibles value to test before the dichotomy
                in order to optimize the computation time.

        Returns:
            float: the approximated root computed.

        """
        hints = [*hints, a, b]
        for x in hints:
            if abs(f(x)) < precision:
                return x

        c = a+ (b-a)/2
        while b-a > precision:
            if abs(f(c)) < precision:
                return c
            if(f(a) > 0 and f(b) < -0):
                if(f(c) >= 0):
                    a = c
                else:
                    b = c
            elif(f(b) > 0 and f(a) < -0):
                if(f(c)>0):
                    b = c
                else:
                    a = c
            else:
                msg = ("Dichotomy is impossible with the function and the domain given"
                       " as parameters. You should check if the function change of sign"
                       " exactly once in the domain.")
                raise ValueError(msg)

            c = a+(b-a)/2

        if abs(f(c)) < abs(f(a)) and abs(f(c)) < abs(f(b)):
            return c
        if abs(f(a)) < abs(f(b)):
            return a
        return b
