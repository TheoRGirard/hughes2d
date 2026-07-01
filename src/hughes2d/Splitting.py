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

import csv
import json
import threading
from collections.abc import Callable
import datetime
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike

from hughes2d.EikonalSolver import EikoSolver
from hughes2d.LWR2D import LWRSolver
from hughes2d.Mesh2D import CellValueMap, Mesh, VertexValueMap

if __name__ == "__main__":
    threading.freeze_support()

previous_thread_dens = object()
previous_thread_vec = object()
EMPTY_LIST = []
DEFAULT_OPTIONS = {"model": "hughes"}
EMPTY_THRESHOLD = 1e-2


class PedestrianSolver:
    """An object wrapping together different models for pedestrian flows.

    Examples:
        First, we need a `Mesh` object in order to create the solver::

            MyMesh = Mesh()
            MyMesh.loadFromJson("pathToMyMesh.json")

        Then we need to declare the options in a dictionnary::

            opt = dict( model = "hughes",
                        filename = "pathForTheSavedFile",
                        save = True,
                        verbose = False,
                        lwrSolver = {   "convexFlux" : True,
                                        "method" : "midVector",
                                        "ApproximationThreshold" : 0.0001},
                        eikoSolver = {  "method" : "FMT",
                                        "constrained" : True,
                                        "NarrowBandDepth" : 2})

        Eventually, we need to define an initial datum and we are in a position to
        instantiate the solver::

            InitialDatum = Mesh2D.CellValueMap(MyMesh)
            InitialDatum.generateRandom()

            MySolver = PedestrianSolver(MyMesh, 0.01, initial_density = InitialDatum, options=opt)

        Now we can compute the approximation until the domain is empty::

            MySolver.compute_until_empty()

    Models:
        At the moment, three models are available for simulations:

          - Hughes" model: a model where the agents try to minimize their individual
            cost by avoiding high density regions.
          - Colombo-Garavello-Lecureux-Mercier-Lecureux-Mercier: a model where the agents try to take
            the shortest path to the exits but are deviated by high density regions.
          - LWR with constant direction field: a model where the agents take the
            shortest path to the exits without taking the surrounding density into
            account.

        See the :doc:`maths` section of the documentation for more details.

    Args:
        Mesh (Mesh): the mesh on which the approximations will be computed.
        dt (float): the duration of a time step. Be careful of the CFL condition
            see :ref:`CFLwarning`.
        initial_density (CellValueMap): the initial density in the domain.
        speed_function (function, float -> float, optional): the speed function
            corresponding to the speed of agents depending on the local density.
        cost_function (function, float -> float, optional): the cost function
            corresponding to the running cost in the eikonal equation. Useless if the
            model used is not "hughes".
        directions (list[list[float]], optional): the direction vector field to use as
            trajectories for the agents. Useless for Hughes" model as the vector field
            is recomputed depending on the density.
            If not prescribed, the vector field is computed at the initialization of
            the solver as the shortest path towards the exits.
        options (dict, optional): an optional dictionary prescribing the model to use
            and various parameters for the numerical simulations.
            See :ref:`options-pedestrian` below.

    .. _options-pedestrian:

    Options:
        Options are passed as an optional dictionary as a parameter of the ``PedestrianSolver`` object.

        ========================== ===== ========================================================== ===========================================================================================================================
        key                        type  possible values                                            description
        ========================== ===== ========================================================== ===========================================================================================================================
        "model"                    str   "hughes", "colombo-garavello" or "constantDirectionField"  determines the model to use for the numerical simulation.
        "save"                     bool  True or False                                              determines whether the data of the simulation should be stored in .csv files.
        "filename"                 str   a valid path                                               sets the path and basename for the save files.
        "framerate"                int   > 0                                                        number of frame per seconds that will be saved in the .csv files. Useless if "save" = False.
        "additional_computations"  dict                                                             adds computations of non standard quantities to the simulation. See :ref:`AdditionnalComputations` for more information.
        "verbose"                  bool  True or False                                              determines if the solver will print informations in the console or not.
        "lwrSolver"                dict                                                             the dictionary containing all the options to use for the ``LWRSolver`` object (see the :ref:`LWRSolver` doc).
        "eikoSolver"               dict                                                             the dictionary containing all the options to use for the ``EikoSolver`` object (see the :ref:`EikoSolver` doc).
        "CGparameters"             dict                                                             the dictionary containing all the parameters to use for the Rinaldo-Garavello-Lecureux-Mercier (see :ref:`CGparameters`)
        ========================== ===== ========================================================== ===========================================================================================================================

        .. _AdditionnalComputations:

        additional_computations
          This dictionary can contain 3 optional keys at the moment:

              - ``total_mass``: computes the total mass in the domain for each time
                step.
              - ``zones_mean_density``: computes the mean density in each zone of the
                mesh defined in ``Mesh.zones`` for each time step.
              - ``max_density``: computes the maximal density present in the domain for
                each time step.

        .. _CGparameters:

        CGparameters
          This dictionary contains two different keys:

              - ``radius`` (float): corresponds to the radius of the convolution.
              - ``epsilon`` (flaoat): corresponds to the amount of influence of the
                deviation vector field.

          See :ref:`ColomboGaravelloModel` in the documention for more details.

    Raises:
        ValueError: if the "model" key in the option dictionary is not set properly.

    Attributes
    ----------

    Attributes:
        mesh (Mesh): the mesh on which the approximations will be computed.
        time_step (int): number of time steps computed since the initialization of the
            solver.
        options (dict): the options dictionary. See :ref:`optionsDict`.
        dt (float): the duration of a time step.
        speed_function (function, float -> float): the speed function corresponding to
            the speed of agents depending on the local density.
        cost_function (function, float -> float): the cost function corresponding to the
            running cost in the eikonal equation. Useless if the model used is not
            "hughes".
        directions (list[list[float]]): the direction vector field to use as
            trajectories for the agents.
        numForgottenSteps (int): number of steps not saved for each time step saved.
            Useless if ``options["save"] == False``.

    Methods
    -------
    """

    def __init__(
        self,
        mesh: Mesh,
        dt: float,
        initial_density: CellValueMap,
        speed_function: Callable[[float], float] = (lambda x: 1 - x),
        cost_function: Callable[[float], float] = (lambda x: 1 + 2 * x),
        directions: list[list[float]] = EMPTY_LIST,
        options: dict = DEFAULT_OPTIONS,
    ) -> None:
        self.mesh = mesh
        self.time_step = 0

        self.options = options

        if "framerate" not in self.options:
            self.options["framerate"] = 25

        self.dt = dt

        if "verbose" not in self.options:
            self.options["verbose"] = False

        self.num_forgotten_steps = max(
            int(1 / (self.dt * self.options["framerate"])), 1
        )
        if self.options["verbose"]:
            print("Number of steps omitted for one frame : ", self.num_forgotten_steps)

        self.speed_function = speed_function
        self.cost_function = cost_function

        if "save" not in self.options:
            self.options["save"] = False

        if "filename" not in self.options:
            self.options["filename"] = "Save" + str(datetime.datetime.now(datetime.UTC))

        if "additional_computations" not in self.options:
            self.options["additional_computations"] = {}

        if "model" not in self.options:
            msg = "A model should be prescribed in the options dictionary."
            raise ValueError(msg)

        if "eikoSolver" not in self.options:
            self.options["eikoSolver"] = {"constrained": True, "NarrowBandDepth": 2}

        if "lwrSolver" not in self.options:
            self.options["lwrSolver"] = {
                "convexFlux": True,
                "anNum": "dichotomy",
                "method": "midVector",
                "ApproximationThreshold": 0.0001,
            }

        if self.options["model"] == "constantDirectionField":
            if len(directions) > 0:
                self.directions = directions
            else:
                self.directions = []
                self.eiko_solver = EikoSolver(
                    self.mesh,
                    density_map=initial_density,
                    cost_function=(lambda x: 1),
                    opt=self.options["eikoSolver"],
                )

                self.eiko_solver.compute_field()

                self.constant_field = self.eiko_solver.field_values

                self.directions = self.eiko_solver.field_values.compute_gradient_flow()

        elif self.options["model"] == "hughes":
            self.directions = []
            self.eiko_solver = EikoSolver(
                self.mesh,
                density_map=initial_density,
                cost_function=self.cost_function,
                opt=self.options["eikoSolver"],
            )

            self.eiko_solver.compute_field()

            self.directions = self.eiko_solver.field_values.compute_gradient_flow()
        elif self.options["model"] == "colombo-garavello":
            if directions != []:
                self.constant_directions = directions
            else:
                self.directions = []
                self.eiko_solver = EikoSolver(
                    self.mesh,
                    density_map=initial_density,
                    cost_function=(lambda x: 1),
                    opt=self.options["eikoSolver"],
                )

                self.eiko_solver.compute_field()

                self.constant_field = self.eiko_solver.field_values

                self.constant_directions = (
                    self.eiko_solver.field_values.compute_gradient_flow()
                )

            self.deviation_field = VertexValueMap(self.mesh)
            if "CGparameters" not in self.options:
                if self.options["verbose"]:
                    print(
                        "Warning: parameters for the Colombo-Garavello-Lecureux-Mercier model not"
                        " prescibed. Defaults as radius = 0.2 and espsilon = 0.1"
                    )
                self.radius = 0.2
                self.epsilon = 0.1
            else:
                if "radius" not in self.options["CGparameters"]:
                    if self.options["verbose"]:
                        print(
                            "Warning: radius for the Colombo-Garavello-Lecureux-Mercier model not"
                            " properly prescibed. Chosing defaults as radius = 0.2"
                        )
                    self.radius = 0.2
                else:
                    self.radius = self.options["CGparameters"]["radius"]

                if "epsilon" not in self.options["CGparameters"]:
                    if self.options["verbose"]:
                        print(
                            "Warning: epsilon for the Colombo-Garavello-Lecureux-Mercier model not"
                            " properly prescibed. Chosing defaults as epsilon = 0.1"
                        )
                    self.epsilon = 0.1
                else:
                    self.epsilon = self.options["CGparameters"]["epsilon"]

            self.normalization_func = lambda x, y: (
                self.radius**2 * (np.sqrt(1 + x**2 + y**2)) / self.epsilon
            )
            self.deviation_field.values = initial_density.convolution_over_square_ball(
                self.radius, self.__convolution_function_cglm
            )

            self.last_density = CellValueMap(self.mesh)

            self.directions = (
                self.constant_directions
                + self.deviation_field.compute_gradient_flow(
                    normalization=self.normalization_func
                )
            )
        else:
            raise ValueError(str(self.options["model"]) + " is not a valid model.")

        if self.options["save"]:
            global previous_thread_dens, previous_thread_vec
            new_thread = threading.Thread(
                target=_write_first_line,
                args=((self.options["filename"] + "_vectors.csv"), self.directions),
            )
            new_thread.start()
            previous_thread_vec = new_thread

            new_thread = threading.Thread(
                target=_write_first_line,
                args=(
                    (self.options["filename"] + "_densities.csv"),
                    initial_density.values,
                ),
            )
            new_thread.start()
            previous_thread_dens = new_thread

        if "total_mass" in self.options["additional_computations"]:
            self.total_masses = [initial_density.integrate()]
        if "zones_mean_density" in self.options["additional_computations"]:
            self.zone_densities = {}
            for zone_name in self.mesh.zones:
                self.zone_densities[zone_name] = []
        if "max_density" in self.options["additional_computations"]:
            self.max_densities = [max(initial_density.values)]

        self.lwr_solver = LWRSolver(
            self.mesh,
            self.dt,
            previous_density=initial_density,
            direction_map=self.directions,
            speed_function=self.speed_function,
            opt=self.options["lwrSolver"],
        )

    def compute_step(self) -> None:
        """Compute one step of time of the approximation of the solution of the chosen model. Also saves the generated data if ``options["save"] == True``."""
        self.time_step += 1
        self.lwr_solver.compute_next_step()
        if self.options["model"] == "constantDirectionField":
            if (self.options["save"]) and (
                self.time_step % self.num_forgotten_steps == 0
            ):
                self.__save_density_slice(self.lwr_solver.densityt1)
        elif self.options["model"] == "hughes":
            self.eiko_solver.update_density(self.lwr_solver.densityt1)
            self.eiko_solver.compute_field()
            self.directions = self.eiko_solver.field_values.compute_gradient_flow()
            if (self.options["save"]) and (
                self.time_step % self.num_forgotten_steps == 0
            ):
                self.__save_density_slice(self.lwr_solver.densityt1)
                self.__save_vector_slice(self.directions)
        elif self.options["model"] == "colombo-garavello":
            self.last_density.values = self.lwr_solver.densityt1
            self.deviation_field.values = (
                self.last_density.convolution_over_square_ball(
                    self.radius, self.__convolution_function_cglm
                )
            )
            self.directions = (
                self.constant_directions
                + self.deviation_field.compute_gradient_flow(
                    normalization=self.normalization_func
                )
            )
            if (self.options["save"]) and (
                self.time_step % self.num_forgotten_steps == 0
            ):
                self.__save_density_slice(self.lwr_solver.densityt1)
                self.__save_vector_slice(self.directions)

        if self.time_step % self.num_forgotten_steps == 0:
            if "total_mass" in self.options["additional_computations"]:
                self.total_masses.append(
                    sum(
                        [
                            self.lwr_solver.densityt1[i] * self.mesh.cell_areas[i]
                            for i in range(len(self.mesh.triangles))
                        ]
                    )
                )
            if "zones_mean_density" in self.options["additional_computations"]:
                for zone_name in self.mesh.zones:
                    self.zone_densities[zone_name].append(
                        sum(
                            [
                                self.lwr_solver.densityt1[i] * self.mesh.cell_areas[i]
                                for i in self.mesh.zones[zone_name]["triangles"]
                            ]
                        )
                        / sum(
                            [
                                self.mesh.cell_areas[i]
                                for i in self.mesh.zones[zone_name]["triangles"]
                            ]
                        )
                    )
            if "max_density" in self.options["additional_computations"]:
                self.max_densities.append(max(self.lwr_solver.densityt1))

        self.lwr_solver.update(self.directions)

    def compute_steps(self, n: int) -> None:
        """Compute ``n`` steps of time of the approximation of the solution of the chosen model.

        Also saves the generated data if ``options["save"] == True``.

        Args:
            n (int): the number of steps to compute.

        """
        if self.options["verbose"]:
            print(f"Computing {n} time steps of the approximated solution.")
        for i in range(n + 1):
            self.compute_step()
            if self.options["verbose"]:
                print("Time step : ", i, "/", n)

        print("End of simulation.")
        self.save_additionnal_computations()

    def compute_until_empty(self, max_frames: int = 5000) -> None:
        """Compute enough steps of time of the approximation so that the domain is empty.

        Also saves the generated data if ``options["save"] == True``.

        Args:
            max_frames (int, optional): the maximal number of steps to compute.

        """
        if self.options["verbose"]:
            print(
                f"Computing time steps until the total mass is less"
                f" than {EMPTY_THRESHOLD}."
                f"Maximum number of steps: {max_frames * self.num_forgotten_steps}."
            )
        if "total_mass" not in self.options["additional_computations"]:
            self.options["additional_computations"]["total_mass"] = True
            self.total_masses = [
                sum(
                    [
                        self.lwr_solver.densityt0[i] * self.mesh.cell_areas[i]
                        for i in range(len(self.mesh.triangles))
                    ]
                )
            ]

        num_steps = 0
        while (
            self.total_masses[-1] > EMPTY_THRESHOLD
            and len(self.total_masses) < max_frames + 1
        ):
            self.compute_step()
            if self.options["verbose"]:
                print(
                    "Time step : ",
                    num_steps,
                    "/",
                    max_frames * self.num_forgotten_steps,
                    " total mass : ",
                    self.total_masses[-1],
                )
            num_steps += 1
        if self.options["verbose"]:
            print("End of simulation.")

        self.save_additionnal_computations()

    def save_additionnal_computations(self) -> None:
        """Save the additional computations required in the ``options["additional_computations"]`` dictionary."""
        if "total_mass" in self.options["additional_computations"]:
            _write_first_line(
                self.options["filename"] + "_total_mass.csv",
                [
                    self.dt * i * self.num_forgotten_steps
                    for i in range(len(self.total_masses))
                ],
            )
            _write_slice(
                self.options["filename"] + "_total_mass.csv", self.total_masses
            )

        if "zones_mean_density" in self.options["additional_computations"]:
            _write_first_line(
                self.options["filename"] + "_zones_mean_density.csv",
                ["time"]
                + [
                    self.dt * i * self.num_forgotten_steps
                    for i in range(
                        len(self.zone_densities[next(iter(self.mesh.zones))])
                    )
                ],
            )
            for zone_name in self.mesh.zones:
                _write_slice(
                    self.options["filename"] + "_zones_mean_density.csv",
                    [zone_name + "_density"] + self.zone_densities[zone_name],
                )

        if "max_density" in self.options["additional_computations"]:
            _write_first_line(
                self.options["filename"] + "_max_density.csv",
                [
                    self.dt * i * self.num_forgotten_steps
                    for i in range(len(self.max_densities))
                ],
            )
            _write_slice(
                self.options["filename"] + "_max_density.csv", self.max_densities
            )

    def compute_steps_and_show(self, n: int) -> None:
        """Compute ``n`` steps of time of the approximation of the solution of the chosen model.

        Then display the solution after ``n`` steps.

        Args:
            n (int): the number of steps to compute.

        """
        if self.options["verbose"]:
            print(f"Computing {n} time steps of the approximated solution.")
        for i in range(n + 1):
            self.compute_step()
            if self.options["verbose"]:
                print("Time step : ", i, "/", n)

        if self.options["model"] != "constantDirectionField":
            self.eiko_solver.field_values.show_vector_field()
        self.lwr_solver.show_density(t=0)

    def __convolution_function_cglm(self, dens: float, dx: float, dy: float) -> float:
        return (
            dens
            * (1 - (dx / self.radius) ** 2) ** 3
            * (1 - (dy / self.radius) ** 2) ** 3
        )

    def __save_density_slice(self, density: ArrayLike) -> None:
        _write_slice_parallel_dens(self.options["filename"] + "_densities.csv", density)

    def __save_vector_slice(self, vector_field: ArrayLike) -> None:
        _write_slice_parallel_vec(
            self.options["filename"] + "_vectors.csv", vector_field
        )

    def save_to_json(self) -> None:
        """Save both the mesh and the informations about the simulation in two separate files '_mesh.json' and '_info.json'."""
        self.mesh.save_to_json(self.options["filename"])

        if self.options["model"] == "constantDirectionField":
            dico = {"type": "density field"}
        else:
            dico = {"type": "vector density field"}
        dico["dt"] = self.dt
        dico["options"] = self.options
        dico["finalTimeStep"] = self.time_step

        dico["densities"] = self.options["filename"] + "_densities.csv"
        dico["potential"] = self.options["filename"] + "_potential.csv"
        dico["vectors"] = self.options["filename"] + "_vectors.csv"

        with Path(self.options["filename"] + "_info.json").open(
            "w", encoding="utf-8"
        ) as f:
            json.dump(dico, f, ensure_ascii=False, indent=4)


def _write_first_line(filename: str, chunk: ArrayLike) -> None:
    with Path(filename).open("w", encoding="UTF8") as f:
        writer = csv.writer(f)
        writer.writerow(chunk)


def _write_slice(filename: str, chunk: ArrayLike) -> None:
    with Path(filename).open("a", encoding="UTF8") as f:
        writer = csv.writer(f)
        writer.writerow(chunk)


def _write_slice_parallel_dens(filename: str, data: ArrayLike) -> None:
    global previous_thread_dens
    if previous_thread_dens.is_alive():
        previous_thread_dens.join()
    new_thread = threading.Thread(target=_write_slice, args=(filename, data))
    new_thread.start()
    previous_thread_dens = new_thread


def _write_slice_parallel_vec(filename: str, data: ArrayLike) -> None:
    global previous_thread_vec
    if previous_thread_vec.is_alive():
        previous_thread_vec.join()

    new_thread = threading.Thread(target=_write_slice, args=(filename, data))
    new_thread.start()
    previous_thread_vec = new_thread
