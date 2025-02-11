
"""
Solver for LWR-like 2D equations

Last update : 29/08/25

Girard Theo
"""


from hughes2d.Mesh2D import *

#import panda as pd
#import plotly
import copy
import json



"""
Objet LWRSolver -uses a transmission map resolution at each edge-
Take as an init parameter :
- A direction map (vector field by triangle)
- A Mesh object
- A density map (cell value map object)
- A speed function (lambda function)
- An boolean optional parameter that tells if the flux function is convex. If this is the case, this parameter makes the computations much lighter.
- An optional parameter prescribing the method of calculation of the approximations:
    dichotomy or ...
- An optional parameter prescribing the numerical method:
    - tmap : transmission map
    - midVector : classical godunov flux with the weighted average vector

Methods :
- computeNextStep : return the computed density for the next time step
- update : set the computed density as a nex initial density, and a new direction field as parameter
"""

class LWRSolver(object):

    def __init__(self, Mesh, dt, dx, previousDensity=[], DirectionMap = [], speedFunction = (lambda x: 1-x), options=dict(convexFlux = True, anNum = "dichotomy", method = "midVector")):
        self.mesh = Mesh #type Mesh
        if(previousDensity != []):
            self.densityt0 = previousDensity #type CellValueMap
        else:
            self.densityt0 = CellValueMap(self.mesh)
        self.densityt1 = CellValueMap(self.mesh)

        #print(self.exitVertices)
        self.dt = dt

        self.dx = dx

        self.options = options

        self.directions = DirectionMap

        if("ApproximationThreshold" not in self.options.keys() ):
            self.options['ApproximationThreshold'] = 0.0001


        self.fluxFunction = lambda x: x*speedFunction(x)

        if(self.options["convexFlux"]):
            if(self.fluxFunction == (lambda x: x*(1-x))):
                self.maxFluxPoint = 0.5
                print("point de max explicite")
            else:
                self.maxFluxPoint = argMax(self.fluxFunction,0,1,self.options['ApproximationThreshold'])
                print("point de max trouvé : ", self.maxFluxPoint)


    def computeStepTmap(self):
        if(self.options["convexFlux"]):
            for triangleIndex, triangleCell in enumerate(self.mesh.trianglesWithEdges):
                Modifdensity = 0
                for edgeNumber, edgeIndex in enumerate(triangleCell):
                    PairOfTriangles = self.mesh.pairsOfTriangles[edgeIndex]
                    TriangleGrad = self.directions[triangleIndex]
                    Normal = self.mesh.outerNormalVectByTriangles[triangleIndex][edgeNumber]

                    triangleFlux = TriangleGrad[0]*Normal[0] + TriangleGrad[1]*Normal[1]

                    if(edgeIndex in self.mesh.exitEdges):
                        farFlux = 1
                        farDensity = 0
                    elif(edgeIndex in self.mesh.wallEdges):
                        triangleFlux = 0
                        farFlux = 0
                        farDensity = 0
                    else:
                        if(PairOfTriangles[0] == triangleIndex):
                            farTriangleGrad = self.directions[PairOfTriangles[1]]
                            farDensity = self.densityt0[PairOfTriangles[1]]
                        else:
                            farTriangleGrad = self.directions[PairOfTriangles[0]]
                            farDensity = self.densityt0[PairOfTriangles[0]]
                        farFlux = (farTriangleGrad[0]*Normal[0] + farTriangleGrad[1]*Normal[1])

                    if(triangleFlux > 1 + self.options['ApproximationThreshold'] or farFlux > 1 + self.options['ApproximationThreshold']):
                        print("ALERT : scalar product out of bonds...")
                        print("Norme farTriangleGrad :", farTriangleGrad[0]*farTriangleGrad[0] + farTriangleGrad[1]*farTriangleGrad[1])
                        print("Norme normal :", Normal[0]*Normal[0] + Normal[1]*Normal[1])


                    """
                    We search for k such that God(triangleFlux function, densityt0[triangleIndex], k ) = God(farFlux function, k, fardensity).
                    In the specific case of a convex flux function we can solve that by treating different cases :
                    """
                    if(triangleFlux <= 0 and farFlux >= 0):
                        Totalflux = 0
                    elif(triangleFlux >= 0 and farFlux <= 0):
                        Totalflux = 0

                    elif(triangleFlux < 0 and farFlux < 0):
                        OutingFlux = triangleFlux*self.fluxFunction(max(self.densityt0[triangleIndex],self.maxFluxPoint))
                        EnteringFlux = farFlux*self.fluxFunction(min(farDensity,self.maxFluxPoint))
                        Totalflux = max(OutingFlux,EnteringFlux)

                    elif(triangleFlux > 0 and farFlux > 0):
                        EnteringFlux = farFlux*self.fluxFunction(max(farDensity,self.maxFluxPoint))
                        OutingFlux = triangleFlux*self.fluxFunction(min(self.densityt0[triangleIndex],self.maxFluxPoint))
                        Totalflux = min(OutingFlux,EnteringFlux)
                    Modifdensity -= (self.dt/self.mesh.cellAreas[triangleIndex]) * self.mesh.edgeLength[edgeIndex]* Totalflux
                self.densityt1[triangleIndex] = min(1,max(self.densityt0[triangleIndex] + Modifdensity,0))

        else: #flux not convex
            for triangleIndex, triangleCell in enumerate(self.mesh.trianglesWithEdges):
                Modifdensity = 0
                for edgeNumber, edgeIndex in enumerate(triangleCell):
                    PairOfTriangles = self.mesh.pairsOfTriangles[edgeIndex]
                    TriangleGrad = self.directions[triangleIndex]
                    Normal = self.mesh.outerNormalVectByTriangles[triangleIndex][edgeNumber]

                    triangleFlux = TriangleGrad[0]*Normal[0] + TriangleGrad[1]*Normal[1]

                    if(edgeIndex in self.mesh.exitEdges):
                        farFlux = 1
                        farDensity = 0
                    elif(edgeIndex in self.mesh.wallEdges):
                        triangleFlux = 0
                        farFlux = 0
                        farDensity = 0
                    else:
                        if(PairOfTriangles[0] == triangleIndex):
                            farTriangleGrad = self.directions[PairOfTriangles[1]]
                            farDensity = self.densityt0[PairOfTriangles[1]]
                        else:
                            farTriangleGrad = self.directions[PairOfTriangles[0]]
                            farDensity = self.densityt0[PairOfTriangles[0]]
                        farFlux = (farTriangleGrad[0]*Normal[0] + farTriangleGrad[1]*Normal[1])

                    if(triangleFlux > 1 + self.options['ApproximationThreshold'] or farFlux > 1 + self.options['ApproximationThreshold']):
                        print("ALERT : scalar product out of bonds...")
                        print("Norme farTriangleGrad :", farTriangleGrad[0]*farTriangleGrad[0] + farTriangleGrad[1]*farTriangleGrad[1])
                        print("Norme normal :", Normal[0]*Normal[0] + Normal[1]*Normal[1])
                        """
                        If the flux function is not convex, we use a dichotomy to solve the problem.
                        """
                    if(triangleFlux <= 0 and farFlux >= 0):
                        Totalflux = 0
                    elif(triangleFlux >= 0 and farFlux <= 0):
                        Totalflux = 0
                    else:
                        parametrizedFlux = lambda k : God(Outerflux, self.densityt0[triangleIndex], k,self.options['ApproximationThreshold']) - God(Innerflux,k, farDensity, self.options['ApproximationThreshold'])
                        k = ApproZeroDichotomie(parametrizedFlux,0,1,self.options['ApproximationThreshold'])
                        Totalflux = God(Outerflux, self.densityt0[triangleIndex], k,self.options['ApproximationThreshold'])
                Modifdensity -= (self.dt/self.mesh.cellAreas[triangleIndex]) * self.mesh.edgeLength[edgeIndex]* Totalflux
            self.densityt1[triangleIndex] = min(1,max(self.densityt0[triangleIndex] + Modifdensity,0))

    def computeStepMidVector(self):
        if(self.options["convexFlux"]):
            for triangleIndex, triangleCell in enumerate(self.mesh.trianglesWithEdges):
                Modifdensity = 0
                for edgeNumber, edgeIndex in enumerate(triangleCell):
                    PairOfTriangles = self.mesh.pairsOfTriangles[edgeIndex]
                    TriangleGrad = self.directions[triangleIndex]
                    Normal = self.mesh.outerNormalVectByTriangles[triangleIndex][edgeNumber]
                    Totalflux = 1

                    if(edgeIndex in self.mesh.exitEdges):
                        VectorFlux = TriangleGrad
                        farDensity = 0
                    elif(edgeIndex in self.mesh.wallEdges):
                        VectorFlux = TriangleGrad
                        Totalflux = 0
                    else:

                        if(PairOfTriangles[0] == triangleIndex):
                            farTriangleGrad = self.directions[PairOfTriangles[1]]
                            farDensity = self.densityt0[PairOfTriangles[1]]
                        else:
                            farTriangleGrad = self.directions[PairOfTriangles[0]]
                            farDensity = self.densityt0[PairOfTriangles[0]]

                        VectorFlux = [(TriangleGrad[0]*self.densityt0[triangleIndex] + farTriangleGrad[0]*farDensity), (TriangleGrad[1]*self.densityt0[triangleIndex] + farTriangleGrad[1]*farDensity) ]
                    normeVectorFlux = np.sqrt(VectorFlux[0]*VectorFlux[0] + VectorFlux[1]*VectorFlux[1])

                    if(normeVectorFlux < self.options['ApproximationThreshold'] or Totalflux == 0):
                        Totalflux = 0
                    else:
                        Scalar = VectorFlux[0]*Normal[0] + VectorFlux[1]*Normal[1]
                        #ModifFluxFunc = lambda x : (VectorFlux[0]*Normal[0] + VectorFlux[1]*Normal[1])/normeVectorFlux *self.fluxFunction(x)
                        #Totalflux = God(ModifFluxFunc, self.densityt0[triangleIndex], farDensity)
                        if(Scalar > 0):
                            if(self.densityt0[triangleIndex] <= farDensity):
                                Totalflux = Scalar*min(self.fluxFunction(self.densityt0[triangleIndex]),self.fluxFunction(farDensity))/normeVectorFlux
                            elif(farDensity < self.maxFluxPoint and self.maxFluxPoint < self.densityt0[triangleIndex]):
                                Totalflux = Scalar*self.fluxFunction(self.maxFluxPoint)/normeVectorFlux
                            else:
                                Totalflux = Scalar*max(self.fluxFunction(self.densityt0[triangleIndex]),self.fluxFunction(farDensity))/normeVectorFlux
                        else:
                            if(self.densityt0[triangleIndex] >= farDensity):
                                Totalflux = Scalar*min(self.fluxFunction(self.densityt0[triangleIndex]),self.fluxFunction(farDensity))/normeVectorFlux
                            elif(farDensity > self.maxFluxPoint and self.maxFluxPoint > self.densityt0[triangleIndex]):
                                Totalflux = Scalar*self.fluxFunction(self.maxFluxPoint)/normeVectorFlux
                            else:
                                Totalflux = Scalar*max(self.fluxFunction(self.densityt0[triangleIndex]),self.fluxFunction(farDensity))/normeVectorFlux

                    Modifdensity -= (self.dt/self.mesh.cellAreas[triangleIndex]) * self.mesh.edgeLength[edgeIndex]* Totalflux
                self.densityt1[triangleIndex] = min(1,max(self.densityt0[triangleIndex] + Modifdensity,0))
            return self.densityt1
        else: #si flux non convexe
            for triangleIndex, triangleCell in enumerate(self.mesh.trianglesWithEdges):
                Modifdensity = 0
                for edgeNumber, edgeIndex in enumerate(triangleCell):
                    PairOfTriangles = self.mesh.pairsOfTriangles[edgeIndex]
                    TriangleGrad = self.directions[triangleIndex]
                    Normal = self.mesh.outerNormalVectByTriangles[triangleIndex][edgeNumber]
                    Totalflux = 1

                    if(edgeIndex in self.mesh.exitEdges):
                        VectorFlux = TriangleGrad
                        farDensity = 0
                    elif(edgeIndex in self.mesh.wallEdges):
                        VectorFlux = TriangleGrad
                        Totalflux = 0
                    else:

                        if(PairOfTriangles[0] == triangleIndex):
                            farTriangleGrad = self.directions[PairOfTriangles[1]]
                            farDensity = self.densityt0[PairOfTriangles[1]]
                        else:
                            farTriangleGrad = self.directions[PairOfTriangles[0]]
                            farDensity = self.densityt0[PairOfTriangles[0]]

                        VectorFlux = [(TriangleGrad[0]*self.densityt0[triangleIndex] + farTriangleGrad[0]*farDensity), (TriangleGrad[1]*self.densityt0[triangleIndex] + farTriangleGrad[1]*farDensity) ]
                    normeVectorFlux = np.sqrt(VectorFlux[0]*VectorFlux[0] + VectorFlux[1]*VectorFlux[1])

                    if(normeVectorFlux < self.options['ApproximationThreshold'] or Totalflux == 0):
                        Totalflux = 0
                    else:
                        ModifFluxFunc = lambda x : (VectorFlux[0]*Normal[0] + VectorFlux[1]*Normal[1])/normeVectorFlux *self.fluxFunction(x)
                        Totalflux = God(ModifFluxFunc, self.densityt0[triangleIndex], farDensity)

                    Modifdensity -= (self.dt/self.mesh.cellAreas[triangleIndex]) * self.mesh.edgeLength[edgeIndex]* Totalflux
                self.densityt1[triangleIndex] = min(1,max(self.densityt0[triangleIndex] + Modifdensity,0))
            return self.densityt1

    def computeNextStep(self):
        if(self.options["method"] == "tmap"):
            return self.computeStepTmap()
        elif(self.options["method"] == "midVector"):
            return self.computeStepMidVector()

    def update(self, newDirectionField):
        self.densityt0 = self.densityt1
        self.directions = newDirectionField

def argMax(f,a,b, SeuilApproMax = 0.0001):
    NumSlice = int(1+1/SeuilApproMax)
    Max = -1000000000
    Pas = abs(b-a)/NumSlice
    xmax = min(a,b)

    for i in range(NumSlice):
        Test = f(min(a,b) + i*Pas)
        if(Test > Max):
            Max = Test
            xmax = min(a,b) + i*Pas
    return(xmax)

def Max(f,a,b, SeuilApproMax = 0.0001):
    Max = -1000000000
    Pas = abs(b-a)/SeuilApproMax

    for i in range(SeuilApproMax):
        Test = f(min(a,b) + i*Pas)
        if(Test > Max):
            Max = Test
    return(Max)

def Min(f,a,b,SeuilApproMax = 0.0001):
    Min = 100000000
    Pas = abs(b-a)/SeuilApproMax

    for i in range(SeuilApproMax):
        Test = f(min(a,b) + i*Pas)
        if(Test < Min):
            Min = Test
    return(Min)

def God(f,a,b, SeuilApproMax = 0.0001):
        if(a < b):
            return Min(f,a,b, SeuilApproMax)
        else :
            return Max(f,b,a, SeuilApproMax)

def ApproZeroDichotomie(f,a,b, SeuilDicho = 0.0001):
    c = a+ (b-a)/2
    while b-a > SeuilDicho:
        #print(a,f(a),b,f(b))
        if(f(a) >= 0 and f(b) <= 0):
            if(f(c) >= 0):
                a = c
            else:
                b = c
        elif(f(b) >= 0 and f(a) <= 0):
            if(f(c)>=0):
                b=c
            else:
                a = c
        else :
            print("impossible de faire une dichotomie sur la fonction donnée")
            print(a,f(a),b,f(b))
            if(abs(f(a)) < SeuilDicho):
                return a
            if(abs(f(b)) < SeuilDicho):
                return b
            return("error")
        c = a+(b-a)/2
        #print(a,c,b)
    return(c)
