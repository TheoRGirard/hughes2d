from hughes2d.Mesh2D import *
from hughes2d.EikonalSolver import *
from hughes2d.LWR2D import *

import multiprocessing

import csv
from datetime import date

previousProcessDens = object()
previousProcessVec = object()

class PedestrianSolver(object):
    """
    An object wrapping together different models for pedestrian flows.

    Examples:
        First, we need a `Mesh` object in order to create the solver::

            MyMesh = Mesh()
            Mesh.loadFromJson("pathToMyMesh.json")

        Then we need to declare the options in a dictionnary::

            opt = dict( model = "hughes",
                        filename = "pathForTheSavedFile",
                        save = True,
                        verbose = False,
                        lwrSolver = {   'convexFlux' : True,
                                        'method' : "midVector",
                                        'ApproximationThreshold' : 0.0001},
                        eikoSolver = {  'method' : "FMT",
                                        'constrained' : True,
                                        'NarrowBandDepth' : 2})

        Eventually, we need to define an initial datum and we are in a position to instantiate the solver::

            InitialDatum = Mesh2D.CellValueMap(MyMesh)
            InitialDatum.generateRandom()

            MySolver = PedestrianSolver(MyMesh, 0.01, initialDensity = InitialDatum, options=opt)

        Now we can compute the approximation until the domain is empty::

            MySolver.computeUntilEmpty()

    Models
    --------

    At the moment, three models are available for simulations:

        - Hughes' model: a model where the agents try to minimize their individual cost by avoiding high density regions.
        - Colombo-Garavello-Lecureux-Mercier: a model where the agents try to take the shortest path to the exits but are deviated by high density regions.
        - LWR with constant direction field: a model where the agents take the shortest path to the exits without taking the surrounding density into account.

    See the :doc:`maths` section of the documentation for more details.

    .. _options-pedestrian:

    Options
    --------

    Options are passed as an optional dictionary as a parameter of the ``PedestrianSolver`` object.

    ========================== ===== ========================================================== ===========================================================================================================================
    key                        type  possible values                                            description
    ========================== ===== ========================================================== ===========================================================================================================================
    'model'                    str   "hughes", "colombo-garavello" or "constantDirectionField"  determines the model to use for the numerical simulation.
    'save'                     bool  True or False                                              determines whether the data of the simulation should be stored in .csv files.
    'filename'                 str   a valid path                                               sets the path and basename for the save files.
    'framerate'                int   > 0                                                        number of frame per seconds that will be saved in the .csv files. Useless if 'save' = False.
    'additional_computations'  dict  -                                                          adds computations of non standard quantities to the simulation. See :ref:`AdditionnalComputations` for more information.
    'verbose'                  bool  True or False                                              determines if the solver will print informations in the console or not.
    'lwrSolver'                dict  -                                                          the dictionary containing all the options to use for the ``LWRSolver`` object (see the :ref:`LWRSolver` doc).
    'eikoSolver'               dict  -                                                          the dictionary containing all the options to use for the ``EikoSolver`` object (see the :ref:`EikoSolver` doc).
    'CGparameters'             dict  -                                                          the dictionary containing all the parameters to use for the Rinaldo-Garavello-Lecureux-Mercier (see :ref:`CGparameters`)
    ========================== ===== ========================================================== ===========================================================================================================================

    .. _AdditionnalComputations:

    additional_computations
        This dictionary can contain 3 optional keys at the moment:

            - ``total_mass``: computes the total mass in the domain for each time step.
            - ``zones_mean_density``: computes the mean density in each zone of the mesh defined in ``Mesh.zones`` for each time step.
            - ``max_density``: computes the maximal density present in the domain for each time step.

    .. _CGparameters:

    CGparameters
        This dictionary contains two different keys:

            - ``radius`` (float): corresponds to the radius of the convolution.
            - ``epsilon`` (flaoat): corresponds to the amount of influence of the deviation vector field.

        See :ref:`ColomboGaravelloModel` in the documention for more details.

    Args:
        Mesh (Mesh): the mesh on which the approximations will be computed.
        dt (float): the duration of a time step. Be careful of the CFL condition see :ref:`CFLwarning`.
        initialDensity (CellValueMap): the initial density in the domain.
        speedFunction (function, float -> float, optional): the speed function corresponding to the speed of agents depending on the local density.
        costFunction (function, float -> float, optional): the cost function corresponding to the running cost in the eikonal equation. Useless if the model used is not "hughes".
        directions (List[List[float]], optional): the direction vector field to use as trajectories for the agents. Useless for Hughes' model as the vector field is recomputed depending on the density.
            If not prescribed, the vector field is computed at the initialization of the solver as the shortest path towards the exits.
        options (dict, optional): an optional dictionary prescribing the model to use and various parameters for the numerical simulations. See :ref:`options-pedestrian` above.

    Raises:
        ValueError: if the "model" key in the option dictionary is not set properly.

    Attributes
    ------------

    Attributes:
        mesh (Mesh): the mesh on which the approximations will be computed.
        timeStep (int): number of time steps computed since the initialization of the solver.
        options (dict): the options dictionary. See :ref:`optionsDict`.
        dt (float): the duration of a time step.
        speedFunction (function, float -> float): the speed function corresponding to the speed of agents depending on the local density.
        costFunction (function, float -> float): the cost function corresponding to the running cost in the eikonal equation. Useless if the model used is not "hughes".
        directions (List[List[float]]): the direction vector field to use as trajectories for the agents.
        numForgottenSteps (int): number of steps not saved for each time step saved. Useless if ``options['save'] == False``.

    Methods
    ----------

    """

    def __init__(self, Mesh:Mesh, dt:float, initialDensity:CellValueMap, speedFunction = (lambda x: 1-x), costFunction = (lambda x: 1+2*x), directions:List[List[float]] = [], options:dict=dict(model="hughes")):
        self.mesh = Mesh
        self.timeStep = 0

        self.options = options

        if('framerate' not in self.options.keys()):
            self.options['framerate'] = 25

        self.dt = dt

        if('verbose' not in self.options.keys()):
            self.options['verbose'] = False

        self.numForgottenSteps = max(int(1/(self.dt*self.options['framerate'])), 1)
        if self.options['verbose']:
            print("Number of steps omitted for one frame : ", self.numForgottenSteps)


        self.speedFunction = speedFunction
        self.costFunction = costFunction

        if('save' not in self.options.keys()):
            self.options['save'] = False

        if('filename' not in self.options.keys()):
            self.options['filename'] = "Save"+str(date.today())

        if('additional_computations' not in self.options.keys()):
            self.options['additional_computations'] = dict()

        if("model" not in self.options.keys()):
            raise ValueError("A model should be prescribed in the options dictionary.")

        if("eikoSolver" not in self.options.keys()):
            self.options['eikoSolver'] = {  'constrained' : True,
                                            'NarrowBandDepth' : 2}

        if("lwrSolver" not in self.options.keys()):
            self.options['lwrSolver'] = {   'convexFlux' : True,
                                            'anNum' : "dichotomy",
                                            'method' : "midVector",
                                            'ApproximationThreshold' : 0.0001}

        if(self.options["model"] == "constantDirectionField"):
            if(len(directions) > 0):
                self.directions = directions
            else:
                self.directions = []
                self.Eikosolver = EikoSolver(self.mesh, DensityMap = initialDensity, costFunction = (lambda x : 1), opt=self.options['eikoSolver'])

                self.Eikosolver.computeField()

                self.constantField = self.Eikosolver.fieldValues

                self.directions = self.Eikosolver.fieldValues.computeGradientFlow()

        elif(self.options["model"] == "hughes"):
            self.directions = []
            self.Eikosolver = EikoSolver(self.mesh, DensityMap = initialDensity, costFunction = self.costFunction, opt=self.options['eikoSolver'])

            self.Eikosolver.computeField()

            self.directions = self.Eikosolver.fieldValues.computeGradientFlow()
        elif(self.options["model"] == "colombo-garavello"):
            if(directions != []):
                self.constantDirections = directions
            else:
                self.directions = []
                self.Eikosolver = EikoSolver(self.mesh, DensityMap = initialDensity, costFunction = (lambda x : 1), opt=self.options['eikoSolver'])

                self.Eikosolver.computeField()

                self.constantField = self.Eikosolver.fieldValues

                self.constantDirections = self.Eikosolver.fieldValues.computeGradientFlow()

            self.deviationField = VertexValueMap(self.mesh)
            if 'CGparameters' not in self.options.keys():
                if self.options['verbose']:
                    print("Warning: parameters for the Colombo-Garavello model not prescibed. Chosing defaults as radius = 0.2 and espsilon = 0.1")
                self.radius = 0.2
                self.epsilon = 0.1
            else:
                if 'radius' not in self.options["CGparameters"].keys():
                    if self.options['verbose']:
                        print("Warning: radius for the Colombo-Garavello model not properly prescibed. Chosing defaults as radius = 0.2")
                    self.radius = 0.2
                else:
                    self.radius = self.options["CGparameters"]["radius"]

                if 'epsilon' not in self.options["CGparameters"].keys():
                    if self.options['verbose']:
                        print("Warning: epsilon for the Colombo-Garavello model not properly prescibed. Chosing defaults as epsilon = 0.1")
                    self.epsilon = 0.1
                else:
                    self.epsilon = self.options["CGparameters"]["epsilon"]

            self.normalizationFunc = (lambda x,y: self.radius**2 * (np.sqrt(1 + x**2 + y**2))/self.epsilon)
            self.deviationField.values = initialDensity.convolutionOverSquareBall(self.radius, self.convolutionFunctionCG)

            self.lastDensity = CellValueMap(self.mesh)

            self.directions = self.constantDirections + self.deviationField.computeGradientFlow(normalization = self.normalizationFunc)
        else:
            raise ValueError( str(self.options["model"]) + " is not a valid model.")

        if(self.options['save']):
            global previousProcessDens, previousProcessVec
            proc = multiprocessing.Process(target=writeFirstLine, args = (self.options['filename']+"_vectors.csv", self.directions))
            proc.start()
            previousProcessVec = proc

            proc = multiprocessing.Process(target=writeFirstLine, args = (self.options['filename']+"_densities.csv", initialDensity.values))
            proc.start()
            previousProcessDens = proc


        if('total_mass' in self.options['additional_computations'].keys()):
            self.totalMass = [initialDensity.integrate()]
        if('zones_mean_density' in self.options['additional_computations'].keys()):
            self.zoneDensity = dict()
            for zoneName in self.mesh.zones.keys():
                self.zoneDensity[zoneName] = []
        if('max_density' in self.options['additional_computations'].keys()):
            self.maxDensity = [max(initialDensity.values)]

        self.LWRsolver = LWRSolver(self.mesh, self.dt, previousDensity = initialDensity, directionMap = self.directions, speedFunction = self.speedFunction, opt = self.options['lwrSolver'])

    def computeStep(self) -> None:
        """
        Computes one step of time of the approximation of the solution of the chosen model. Also saves the generated data if ``options['save'] == True``.
        """
        self.timeStep += 1
        self.LWRsolver.computeNextStep()
        if(self.options["model"] == "constantDirectionField"):
            if(self.options['save']):
                if(self.timeStep % self.numForgottenSteps == 0):
                    self.saveDensityslice(self.LWRsolver.densityt1)
        elif(self.options["model"] == "hughes"):
            self.Eikosolver.updateDensity(self.LWRsolver.densityt1)
            self.Eikosolver.computeField()
            self.directions = self.Eikosolver.fieldValues.computeGradientFlow()
            if(self.options['save']):
                if(self.timeStep % self.numForgottenSteps == 0):
                    self.saveDensityslice(self.LWRsolver.densityt1)
                    self.saveVectorslice(self.directions)
        elif(self.options["model"] == "colombo-garavello"):
            self.lastDensity.values = self.LWRsolver.densityt1
            self.deviationField.values = self.lastDensity.convolutionOverSquareBall(self.radius, self.convolutionFunctionCG)
            self.directions = self.constantDirections + self.deviationField.computeGradientFlow(normalization = self.normalizationFunc)
            if(self.options['save']):
                if(self.timeStep % self.numForgottenSteps == 0):
                    self.saveDensityslice(self.LWRsolver.densityt1)
                    self.saveVectorslice(self.directions)

        if(self.timeStep % self.numForgottenSteps == 0):
            if('total_mass' in self.options['additional_computations'].keys()):
                self.totalMass.append(sum([self.LWRsolver.densityt1[i]*self.mesh.cellAreas[i] for i in range(len(self.mesh.triangles))]))
            if('zones_mean_density' in self.options['additional_computations'].keys()):
                for zoneName in self.mesh.zones.keys():
                    self.zoneDensity[zoneName].append(sum([self.LWRsolver.densityt1[i]*self.mesh.cellAreas[i] for i in self.mesh.zones[zoneName]['triangles']])/sum([self.mesh.cellAreas[i] for i in self.mesh.zones[zoneName]['triangles']]))
            if('max_density' in self.options['additional_computations'].keys()):
                self.maxDensity.append(max(self.LWRsolver.densityt1))

        self.LWRsolver.update(self.directions)

    def computeSteps(self, n:int) -> None:
        """
        Computes ``n`` steps of time of the approximation of the solution of the chosen model. Also saves the generated data if ``options['save'] == True``.

        Args:
            n (int): the number of steps to compute.
        """
        for i in range(n):
            self.computeStep()
            if(self.options['verbose']):
                print("Time step : ", i, "/", n)

        self.saveAdditionnalComputations()

    def computeUntilEmpty(self, max_frames = 5000) -> None:
        """
        Computes enough steps of time of the approximation so that the domain is empty. Also saves the generated data if ``options['save'] == True``.

        Args:
            max_frames (int, optional): the maximal number of steps to compute.
        """
        if('total_mass' not in self.options['additional_computations'].keys()):
            self.options['additional_computations']['total_mass'] = True
            self.totalMass = [sum([self.LWRsolver.densityt0[i]*self.mesh.cellAreas[i] for i in range(len(self.mesh.triangles))])]

        numSteps = 0
        while(self.totalMass[-1] > float(1e-2) and len(self.totalMass) < max_frames):
            self.computeStep()
            if(self.options['verbose']):
                print("Time step : ", numSteps, "/", max_frames*self.numForgottenSteps," total mass : ", self.totalMass[-1])
            numSteps += 1

        self.saveAdditionnalComputations()

    def saveAdditionnalComputations(self) -> None:
        """
        Saves the additional computations required in the ``options['additional_computations']`` dictionary.
        """
        if('total_mass' in self.options['additional_computations'].keys()):
            writeFirstLine(self.options['filename']+"_total_mass.csv",[self.dt*i*self.numForgottenSteps for i in range(len(self.totalMass))])
            writeSlice(self.options['filename']+"_total_mass.csv",self.totalMass)

        if('zones_mean_density' in self.options['additional_computations'].keys()):
            writeFirstLine(self.options['filename']+"_zones_mean_density.csv",["time"]+[self.dt*i*self.numForgottenSteps for i in range(len(self.zoneDensity[list(self.mesh.zones.keys())[0]]))])
            for zoneName in self.mesh.zones.keys():
                writeSlice(self.options['filename']+"_zones_mean_density.csv",[zoneName+"_density"]+ self.zoneDensity[zoneName])

        if('max_density' in self.options['additional_computations'].keys()):
            writeFirstLine(self.options['filename']+"_max_density.csv",[self.dt*i*self.numForgottenSteps for i in range(len(self.maxDensity))])
            writeSlice(self.options['filename']+"_max_density.csv",self.maxDensity)

    def computeStepsAndShow(self,n:int) -> None:
        """
        Computes ``n`` steps of time of the approximation of the solution of the chosen model. Then display the solution after ``n`` steps.

        Args:
            n (int): the number of steps to compute.
        """
        for i in range(n):
            self.computeStep()
            if(self.options['verbose']):
                print("Time step : ", i, "/", n)

        if(not self.options['model'] == "constantDirectionField"):
            self.Eikosolver.fieldValues.showVectorField()
        self.LWRsolver.showDensity(t=0)

    def convolutionFunctionCG(self, dens:float, dx:float, dy:float) -> float:
        return dens*(1-(dx/self.radius)**2)**3 * (1-(dy/self.radius)**2)**3

    def saveDensityslice(self, density:ArrayLike) -> None:
        writeSlice_parallel_Dens(self.options['filename']+ "_densities.csv", self.timeStep, density)

    def saveVectorslice(self, vectorField:ArrayLike) -> None:
        writeSlice_parallel_Vec(self.options['filename']+ "_vectors.csv", self.timeStep, vectorField)

    def saveToJson(self):
        """
        Saves both the mesh and the informations about the simulation in two separate files "_mesh.json" and "_info.json"
        """
        self.mesh.saveToJson(self.options['filename'])

        if(self.options["model"] == "constantDirectionField"):
            dico = {"type":"density field"}
        else:
            dico = {"type":"vector density field"}
        dico["dt"] = self.dt
        dico["options"] = self.options
        dico["finalTimeStep"] = self.timeStep

        DensityFilename = self.options['filename'] + "_densities.csv"
        dico["densities"] = self.options['filename'] + "_densities.csv"

        PotentialFilename = self.options['filename'] + "_potential.csv"
        VectorsFilename = self.options['filename'] + "_vectors.csv"
        dico["potential"] = PotentialFilename
        dico["vectors"] = VectorsFilename

        with open(self.options['filename']+"_info.json", 'w', encoding='utf-8') as f:
            json.dump(dico, f, ensure_ascii=False, indent=4)



def writeFirstLine(filename, chunk):
    with open(filename, 'w', encoding='UTF8') as f:
        writer = csv.writer(f)
        writer.writerow(chunk)

def writeSlice(filename, chunk):
    with open(filename, 'a', encoding='UTF8') as f:
        writer = csv.writer(f)
        writer.writerow(chunk)

def writeSlice_parallel_Dens(filename, numSlice, data, num_processes=4):
    global previousProcessDens
    if(previousProcessDens.is_alive()):
        previousProcessDens.join()
    proc = multiprocessing.Process(target=writeSlice, args = (filename,data))
    proc.start()
    previousProcessDens = proc


def writeSlice_parallel_Vec(filename, numSlice, data, num_processes=4):
    global previousProcessVec
    if(previousProcessVec.is_alive()):
        previousProcessVec.join()

    proc = multiprocessing.Process(target=writeSlice, args = (filename,data))
    proc.start()
    previousProcessVec = proc
