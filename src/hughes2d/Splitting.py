from hughes2d.Mesh2D import *
from hughes2d.EikonalSolver import *
from hughes2d.LWR2D import *

import multiprocessing

import csv
from datetime import date

previousProcessDens = object()
previousProcessVec = object()

"""
LWRscheme
A complete solver with multiple options:
- If constantDirectionField is set to True, the user must provide a vector field by triangle as directions parameter.
- Else, the direction field is computed by an Eikonal equation at each time step using the costFunction.
The options dictionnary is passed to the LWRsolver.

Methods:
- ComputeSteps : n is the number of steps to compute
- computeStepsAndShow : compute n step and plot the last time step
- saveToJson : save the results in filename.json
- show : plot the densities
"""

class PedestrianSolver(object):

    def __init__(self, Mesh, dt, dx, initialDensity=[], speedFunction = (lambda x: 1-x), costFunction = (lambda x: 1+2*x), directions = [], options=dict(constantDirectionField = True, convexFlux = True, anNum = "dichotomy", method = "midVector")):
        self.mesh = Mesh #type Mesh
        self.timeStep = 0

        self.options = options

        if('framerate' not in self.options.keys()):
            self.options['framerate'] = 25

        self.dt = dt

        self.numForgottenSteps = max(int(1/(self.dt*self.options['framerate'])), 1)
        print("Number of steps omitted for one frame : ", self.numForgottenSteps)

        self.dx = dx

        self.speedFunction = speedFunction
        print(self.speedFunction(0))
        self.costFunction = costFunction

        if('save' not in self.options.keys()):
            self.options['save'] = False

        if('filename' not in self.options.keys()):
            self.options['filename'] = "Save"+str(date.today())

        if('additional_computations' not in self.options.keys()):
            self.options['additional_computations'] = dict()

        if("model" not in self.options.keys()):
            raise ValueError("A model should be prescribed")

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
            self.radius = self.options["CGparameters"]["radius"]
            self.epsilon = self.options["CGparameters"]["epsilon"]

            self.normalizationFunc = (lambda x,y: self.radius**2 * (np.sqrt(1 + x**2 + y**2))/self.epsilon)
            self.deviationField.values = initialDensity.convolutionOverSquareBall(self.radius, self.convolutionFunctionCG)

            self.lastDensity = CellValueMap(self.mesh)

            self.directions = self.constantDirections + self.deviationField.computeGradientFlow(normalization = self.normalizationFunc)

        if(self.options['save']):
            global previousProcessDens, previousProcessVec
            proc = multiprocessing.Process(target=writeFirstLine, args = (self.options['filename']+"_vectors.csv", self.directions))
            print(self.directions)
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

        self.LWRsolver = LWRSolver(self.mesh, self.dt, self.dx, previousDensity = initialDensity, DirectionMap = self.directions, speedFunction = self.speedFunction, options = self.options['lwrSolver'])

    def computeStep(self):
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

    def computeSteps(self, n):
        for i in range(n):
            self.computeStep()
            if("verbose" in self.options.keys()):
                print("Time step : ", i, "/", n)

        self.saveAdditionnalComputations()

    def computeUntilEmpty(self, max_frames = 5000):
        if('total_mass' not in self.options['additional_computations'].keys()):
            self.options['additional_computations']['total_mass'] = True
            self.totalMass = [sum([self.LWRsolver.densityt0[i]*self.mesh.cellAreas[i] for i in range(len(self.mesh.triangles))])]
            print("première masse trouvée : ", self.totalMass[0])


        numSteps = 0
        while(self.totalMass[-1] > float(1e-2) and len(self.totalMass) < max_frames):
            self.computeStep()
            if("verbose" in self.options.keys()):
                print("Time step : ", numSteps, "/", max_frames*self.numForgottenSteps," total mass : ", self.totalMass[-1])
            numSteps += 1

        self.saveAdditionnalComputations()

    def saveAdditionnalComputations(self):
        if('total_mass' in self.options['additional_computations'].keys()):
            writeFirstLine(self.options['filename']+"_total_mass.csv",[self.dt*i*self.numForgottenSteps for i in range(len(self.totalMass))])
            writeSlice(self.options['filename']+"_total_mass.csv",self.totalMass)

        if('zones_mean_density' in self.options['additional_computations'].keys()):
            writeFirstLine(self.options['filename']+"_zones_mean_densiy.csv",["time"]+[self.dt*i*self.numForgottenSteps for i in range(len(self.zoneDensity[list(self.mesh.zones.keys())[0]]))])
            for zoneName in self.mesh.zones.keys():
                writeSlice(self.options['filename']+"_zones_mean_densiy.csv",[zoneName+"_density"]+ self.zoneDensity[zoneName])

        if('max_density' in self.options['additional_computations'].keys()):
            writeFirstLine(self.options['filename']+"_max_density.csv",[self.dt*i*self.numForgottenSteps for i in range(len(self.maxDensity))])
            writeSlice(self.options['filename']+"_max_density.csv",self.maxDensity)

    def computeStepsAndShow(self,n):
        for i in range(n):
            self.computeStep()
            if("verbose" in self.options.keys()):
                print("Time step : ", i, "/", n)

        if(not self.options["constantDirectionField"]):
            self.directions.showVectorField()
        self.LWRsolver.densityt1.show()

    def convolutionFunctionCG(self, dens, dx, dy):
        return dens*(1-(dx/self.radius)**2)**3 * (1-(dy/self.radius)**2)**3


    def saveDensityslice(self, density):
        writeSlice_parallel_Dens(self.options['filename']+ "_densities.csv", self.timeStep, density)

    def saveVectorslice(self, vectorField):
        writeSlice_parallel_Vec(self.options['filename']+ "_vectors.csv", self.timeStep, vectorField)

    def saveToJson(self):
        if(self.options["constantDirectionField"]):
            dico = {"type":"density field"}
        else:
            dico = {"type":"vector density field"}
        dico["dt"] = self.dt
        dico["options"] = self.options
        dico["finalTimeStep"] = self.timeStep
        self.mesh.appendDict(dico)

        DensityFilename = self.options['filename'] + "_densities.csv"
        dico["densities"] = self.options['filename'] + "_densities.csv"

        PotentialFilename = self.options['filename'] + "_potential.csv"
        VectorsFilename = self.options['filename'] + "_vectors.csv"
        dico["potential"] = PotentialFilename
        dico["vectors"] = VectorsFilename

        with open(self.options['filename']+".json", 'w', encoding='utf-8') as f:
            json.dump(dico, f, ensure_ascii=False, indent=4)



def writeFirstLine(filename, chunk):
    with open(filename, 'w', encoding='UTF8') as f:
        # create the csv writer
        writer = csv.writer(f)
        # write a row to the csv file
        writer.writerow(chunk)

def writeSlice(filename, chunk):
    # open the file in the write mode
    with open(filename, 'a', encoding='UTF8') as f:
        # create the csv writer
        writer = csv.writer(f)
        # write a row to the csv file
        writer.writerow(chunk)

def writeSlice_parallel_Dens(filename, numSlice, data, num_processes=4):
    global previousProcessDens
    if(previousProcessDens.is_alive()):
        previousProcessDens.join()
    proc = multiprocessing.Process(target=writeSlice, args = (filename,data))
    proc.start()
    #proc.join()
    previousProcessDens = proc
    #print("Enregistrement dans ", filename, " de l'étape ", numSlice)

def writeSlice_parallel_Vec(filename, numSlice, data, num_processes=4):
    global previousProcessVec
    if(previousProcessVec.is_alive()):
        previousProcessVec.join()

    proc = multiprocessing.Process(target=writeSlice, args = (filename,data))
    proc.start()
    #proc.join()
    previousProcessVec = proc

"""def show(self):


    fig = go.Figure()
    self.mesh.domain.addPlot(fig)

    #import plotly.graph_objects as go



    for j,T in enumerate(self.mesh.triangles):
        fig.add_trace(go.Scatter(x=[self.mesh.vertices[i][0] for i in T]+[self.mesh.vertices[T[0]][0]],
                                y=[self.mesh.vertices[i][1] for i in T]+[self.mesh.vertices[T[0]][1]],
                        fill="toself",
                        hoverinfo = "none",
                        showlegend = False,
                        mode="none",
                        fillcolor ='rgb('+str( int(255*min(1,max(self.densities[0].values[j],0))) )+',0,0)'
                        ))


    frames = []
    frames.append({'data':copy.deepcopy(fig['data']),'name':f'frame{0}'})



    for dens in self.densities:
        triangleTraces  = []
        for j,T in enumerate(self.mesh.triangles):
            triangleTraces.append(go.Scatter(x=[self.mesh.vertices[i][0] for i in T]+[self.mesh.vertices[T[0]][0]],
                                    y=[self.mesh.vertices[i][1] for i in T]+[self.mesh.vertices[T[0]][1]],
                            fill="toself",
                            hoverinfo = "none",
                            showlegend = False,
                            mode="none",
                            fillcolor ='rgb('+str( int(255*min(1,max(dens.values[j],0))) )+',0,0)'
                            ))
        frames.append(go.Frame(data=triangleTraces))

    fig.update(frames=frames)

    n_frames = len(self.densities)

    updatemenus = [dict(
        buttons = [
                dict(
                    args = [None, {"frame": {"duration": 500, "redraw": True},
                                    "fromcurrent": True}],
                    label = "Play",
                    method = "animate"
                    ),
                dict(
                     args = [[None], {"frame": {"duration": 0, "redraw": False},
                                      "mode": "immediate",
                                      "transition": {"duration": 0}}],
                    label = "Pause",
                    method = "animate"
                    )
            ],
            direction = "left",
            pad = {"r": 10, "t": 87},
            showactive = False,
            type = "buttons",
            x = 0.1,
            xanchor = "right",
            y = 0,
            yanchor = "top"
    )]

    sliders = [dict(steps = [dict(method= 'animate',
                                  args= [[f'frame{k}'],
                                  dict(mode= 'immediate',
                                       frame= dict(duration=400, redraw=True),
                                       transition=dict(duration= 0))
                                     ],
                                  label=f'{k+1}'
                                 ) for k in range(n_frames)],
                    active=0,
                    transition= dict(duration= 0 ),
                    x=0, # slider starting position
                    y=0,
                    currentvalue=dict(font=dict(size=12),
                                      prefix='frame: ',
                                      visible=True,
                                      xanchor= 'center'
                                     ),
                    len=1.0) #slider length
               ]
    fig.update_layout(width=600, height=600,

                      updatemenus=updatemenus,
                      sliders=sliders)

    fig.show()"""
