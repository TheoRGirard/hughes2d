#from Mesh2D import *
from hughes2d.Mesh2D import *
#import plotly
import random as alea


class EikoSolver(object):
    """
    An object wrapping together all the methods used to solve the eikonal equation in a discontinuous context.
    We use the following notations for the eikonal equation:

    .. math::
        |\\nabla u | = c(\\rho).


    Attributes
    ------------

    Attributes:
        mesh (Mesh): the mesh on which the solution to the eikonal equation will be approximated.
        density (CellValueMap): a function constant on the triangle of the mesh representing the density of pedestrians.
        fieldValues (VertexValueMap): the object corresponding to the approximated solution of the eikonal equation.
        cost (function, float -> float): the function *c* that must be applied to the density map. If we set :math:`c: \\rho \\mapsto 1` the solution of the eikonal equation will correspond to the length of the shortest path to the exits, without taking the density into account.
        opt (dict): the dictionary containing all the parameters of computation:

            ================= ====== ================ ===================================================================================================================
            key               type   possible values  description
            ================= ====== ================ ===================================================================================================================
            'method'           str    "FME" or "FMT"   corresponds to the computational method to use
            'NarrowBandDepth'  int    1 or 2           the maximal degree of the neighbours explored at each iteration (only relevant if method="FMT")
            'constrained'      bool   True or False    determines if the gradient computed must be constrained inside the triangle or not (only relevant if method="FMT")
            'debugging'        bool   True or False    debugging option, for verbose outputs and detailled prints
            ================= ====== ================ ===================================================================================================================

    Note:
        The best combination of options with respect to the quality of the approximation appears to be::

            {
                'method' : "FMT",
                'NarrowBandDepth' : 2,
                'constrained' : True
            }

        For further details about the algorithms and their comparison, we defer to https://hal.science/tel-05275359v1

    Methods
    ---------
    """

    def __init__(self, Mesh: Mesh, DensityMap: CellValueMap = [], costFunction = (lambda x: 1+x), opt:dict = dict()):
        self.mesh = Mesh #type Mesh
        if(DensityMap != []):
            self.density = DensityMap #type CellValueMap
        else:
            self.density = CellValueMap(self.mesh)

        #self.fieldValues = np.array([0.0 for _ in self.mesh.vertices])
        self.fieldValues = VertexValueMap(self.mesh)

        self.opt = opt
        if('method' not in self.opt.keys()):
            self.opt['method'] = "FMT"

        if('NarrowBandDepth' not in self.opt.keys()):
            self.opt['NarrowBandDepth'] = 2
        if('constrained' not in self.opt.keys()):
            self.opt['constrained'] = True

        if("debugging" not in self.opt.keys()):
            self.opt['debugging'] = False

        self.cost = costFunction

    def updateDensity(self, density: CellValueMap):
        """
        Updates the density of the solver with the new density passed as a parameter.

        Args:
            density (CellValueMap): the new density to update in the solver.
        """
        self.density = density

    def computeFieldUnconstrainedDep2(self) -> None:
        """
        Computes the approximated solution of the eikonal equation and stores the approximation in 'fieldValues'.
        The method used is the FMT algorithm with no constraint on the gradient and vertices at distance one or two are considered neighbour.
        i.e. constrained = False and NarrowBandDepth = 2 see 'opt'
        """
        StateMap = [0 for i in range(len(self.mesh.vertices))]

        VisitedDict = dict()
        Visited = []
        TriangleOfVisit = []
        Validated = []
        self.fieldValues.setInfinity()


        #Dirichlet conditions on the exits
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            StateMap[index] = 1


        #Construction of the narrow band + computation of the potential values
        for index in Validated:
            for [triangleindex, triangleWithoutIndex] in self.mesh.trianglesPerVertex[index]:

                if(StateMap[triangleWithoutIndex[1]] == 0):
                    StateMap[triangleWithoutIndex[1]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 0):
                    StateMap[triangleWithoutIndex[0]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 1 and StateMap[triangleWithoutIndex[1]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradUnconstrained(self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[index],
                                                                            self.fieldValues[triangleWithoutIndex[0]],
                                                                            self.fieldValues[index],
                                                                            self.cost(self.density[triangleindex]))

                    if(triangleWithoutIndex[1] in VisitedDict.keys()):
                        if(triangleindex in VisitedDict[triangleWithoutIndex[1]]):
                            if VisitedDict[triangleWithoutIndex[1]][triangleindex] == self.fieldValues[triangleWithoutIndex[1]]:
                                self.fieldValues[triangleWithoutIndex[1]] = min([ VisitedDict[triangleWithoutIndex[1]][tridex] for tridex in VisitedDict[triangleWithoutIndex[1]].keys() if tridex != triangleindex]+[PotentialValue])
                                Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[1])
                        VisitedDict[triangleWithoutIndex[1]][triangleindex] = PotentialValue
                    else:
                        VisitedDict[triangleWithoutIndex[1]] = dict()
                        VisitedDict[triangleWithoutIndex[1]][triangleindex] = PotentialValue


                    if PotentialValue < self.fieldValues[triangleWithoutIndex[1]]:
                        self.fieldValues[triangleWithoutIndex[1]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[1])


                elif(StateMap[triangleWithoutIndex[1]] == 1 and StateMap[triangleWithoutIndex[0]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradUnconstrained(self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[index],
                                                                            self.fieldValues[triangleWithoutIndex[1]],
                                                                            self.fieldValues[index],
                                                                            self.cost(self.density[triangleindex]))
                    if(triangleWithoutIndex[0] in VisitedDict.keys()):
                        if(triangleindex in VisitedDict[triangleWithoutIndex[0]]):
                            if VisitedDict[triangleWithoutIndex[0]][triangleindex] == self.fieldValues[triangleWithoutIndex[0]]:
                                self.fieldValues[triangleWithoutIndex[0]] = min([ VisitedDict[triangleWithoutIndex[0]][tridex] for tridex in VisitedDict[triangleWithoutIndex[0]].keys() if tridex != triangleindex]+[PotentialValue])
                                Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[0])
                        VisitedDict[triangleWithoutIndex[0]][triangleindex] = PotentialValue
                    else:
                        VisitedDict[triangleWithoutIndex[0]] = dict()
                        VisitedDict[triangleWithoutIndex[0]][triangleindex] = PotentialValue

                    if PotentialValue < self.fieldValues[triangleWithoutIndex[0]]:
                        self.fieldValues[triangleWithoutIndex[0]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[0])

                elif(StateMap[triangleWithoutIndex[1]] == 2 and StateMap[triangleWithoutIndex[0]] == 2):

                    PotentialValue = self.fieldValues[index] + self.cost(self.density[triangleindex])*EikoSolver.computeHeightLength(self.mesh.vertices[index],
                                                                                                    self.mesh.vertices[triangleWithoutIndex[0]],
                                                                                                    self.mesh.vertices[triangleWithoutIndex[1]])

                    for i in range(2):
                        if( triangleWithoutIndex[i] not in VisitedDict.keys()):
                            VisitedDict[triangleWithoutIndex[i]] = dict()

                        VisitedDict[triangleWithoutIndex[i]][triangleindex] = PotentialValue
                        if PotentialValue < self.fieldValues[triangleWithoutIndex[i]]:
                            #VisitedDict[triangleWithoutIndex[i]][triangleindex]['prev'] = self.fieldValues[triangleWithoutIndex[i]] # [ VertexIndex, triangleIndex,previousBest value ]
                            self.fieldValues[triangleWithoutIndex[i]] = PotentialValue
                            Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[i])


        #Main loop : validation of a vertex and recomputation of the narrow band
        LastValidated = Validated[0]
        NumVertices = len(StateMap)

        while(len(Validated) < NumVertices):

            #Validation of a vertex
            if(len(Visited) == 0):
                if(self.opt['debugging']):
                    print("Completion : ", len(Validated), "/", NumVertices)
                raise ValueError("The fast marching ended too fast... Is the Mesh connected ?")
            if(self.fieldValues[Visited[0]] == float('inf')):
                if(self.opt['debugging']):
                    print("Infinity value detected at vertex (",self.mesh.vertices[Visited[0]][0],",",self.mesh.vertices[Visited[0]][1], ")" )
                raise ValueError("Gradient computation impossible.")

            Validated.append(Visited[0])
            LastValidated = Visited[0]
            Visited = Visited[1:]
            StateMap[LastValidated] = 1

            #recomputation of the narrow band
            for [triangleindex, triangleWithoutIndex] in self.mesh.trianglesPerVertex[LastValidated]:

                if(StateMap[triangleWithoutIndex[1]] == 0):
                    StateMap[triangleWithoutIndex[1]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 0):
                    StateMap[triangleWithoutIndex[0]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 1 and StateMap[triangleWithoutIndex[1]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradUnconstrained(self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[LastValidated],
                                                                            self.fieldValues[triangleWithoutIndex[0]],
                                                                            self.fieldValues[LastValidated],
                                                                            self.cost(self.density[triangleindex]))
                    if(triangleWithoutIndex[1] in VisitedDict.keys()): #Update visited dict
                        if(triangleindex in VisitedDict[triangleWithoutIndex[1]]):
                            if VisitedDict[triangleWithoutIndex[1]][triangleindex] == self.fieldValues[triangleWithoutIndex[1]]:
                                self.fieldValues[triangleWithoutIndex[1]] = min([ VisitedDict[triangleWithoutIndex[1]][tridex] for tridex in VisitedDict[triangleWithoutIndex[1]].keys() if tridex != triangleindex]+[PotentialValue])
                                Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[1])
                        VisitedDict[triangleWithoutIndex[1]][triangleindex] = PotentialValue
                    else:
                        VisitedDict[triangleWithoutIndex[1]] = dict()
                        VisitedDict[triangleWithoutIndex[1]][triangleindex] = PotentialValue

                    if PotentialValue < self.fieldValues[triangleWithoutIndex[1]]:
                        self.fieldValues[triangleWithoutIndex[1]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited, triangleWithoutIndex[1])


                elif(StateMap[triangleWithoutIndex[1]] == 1 and StateMap[triangleWithoutIndex[0]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradUnconstrained(self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[LastValidated],
                                                                            self.fieldValues[triangleWithoutIndex[1]],
                                                                            self.fieldValues[LastValidated],
                                                                            self.cost(self.density[triangleindex]))
                    if(triangleWithoutIndex[0] in VisitedDict.keys()): #Update visited dict
                        if(triangleindex in VisitedDict[triangleWithoutIndex[0]]):
                            if VisitedDict[triangleWithoutIndex[0]][triangleindex] == self.fieldValues[triangleWithoutIndex[0]]:
                                self.fieldValues[triangleWithoutIndex[0]] = min([ VisitedDict[triangleWithoutIndex[0]][tridex] for tridex in VisitedDict[triangleWithoutIndex[0]].keys() if tridex != triangleindex]+[PotentialValue])
                                Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[0])
                        VisitedDict[triangleWithoutIndex[0]][triangleindex] = PotentialValue
                    else:
                        VisitedDict[triangleWithoutIndex[0]] = dict()
                        VisitedDict[triangleWithoutIndex[0]][triangleindex] = PotentialValue

                    if PotentialValue < self.fieldValues[triangleWithoutIndex[0]]:
                        self.fieldValues[triangleWithoutIndex[0]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited, triangleWithoutIndex[0])

                elif(StateMap[triangleWithoutIndex[1]] == 2 and StateMap[triangleWithoutIndex[0]] == 2):

                    PotentialValue = self.fieldValues[LastValidated] + self.cost(self.density[triangleindex])*EikoSolver.computeHeightLength(self.mesh.vertices[LastValidated],
                                                                                                    self.mesh.vertices[triangleWithoutIndex[0]],
                                                                                                    self.mesh.vertices[triangleWithoutIndex[1]])


                    for i in range(2):
                        if( triangleWithoutIndex[i] not in VisitedDict.keys()):
                            VisitedDict[triangleWithoutIndex[i]] = dict()

                        VisitedDict[triangleWithoutIndex[i]][triangleindex] = PotentialValue
                        if PotentialValue < self.fieldValues[triangleWithoutIndex[i]]:
                            self.fieldValues[triangleWithoutIndex[i]] = PotentialValue
                            Visited = self.addInOrderByFieldValue(Visited, triangleWithoutIndex[i])

    def computeFieldConstrainedDep2(self) -> None:
        """
        Computes the approximated solution of the eikonal equation and stores the approximation in 'fieldValues'.
        The method used is the FMT algorithm with constraints on the gradient and vertices at distance one or two are considered neighbour.
        i.e. constrained = True and NarrowBandDepth = 2 see 'opt'
        """
        StateMap = [0 for i in range(len(self.mesh.vertices))]

        Visited = []
        Validated = []
        self.fieldValues.setInfinity()

        #Dirichlet free exit conditions
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            StateMap[index] = 1

        #Construction of the narrow band
        for index in Validated:
            for [triangleindex, triangleWithoutIndex] in self.mesh.trianglesPerVertex[index]:

                if(StateMap[triangleWithoutIndex[1]] == 0):
                    StateMap[triangleWithoutIndex[1]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 0):
                    StateMap[triangleWithoutIndex[0]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 1 and StateMap[triangleWithoutIndex[1]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradConstrained(self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[index],
                                                                            self.fieldValues[triangleWithoutIndex[0]],
                                                                            self.fieldValues[index],
                                                                            self.cost(self.density[triangleindex]))
                    if PotentialValue < self.fieldValues[triangleWithoutIndex[1]]:
                        self.fieldValues[triangleWithoutIndex[1]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[1])


                elif(StateMap[triangleWithoutIndex[1]] == 1 and StateMap[triangleWithoutIndex[0]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradConstrained(self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[index],
                                                                            self.fieldValues[triangleWithoutIndex[1]],
                                                                            self.fieldValues[index],
                                                                            self.cost(self.density[triangleindex]))
                    if PotentialValue < self.fieldValues[triangleWithoutIndex[0]]:
                        self.fieldValues[triangleWithoutIndex[0]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[0])

                elif(StateMap[triangleWithoutIndex[1]] == 2 and StateMap[triangleWithoutIndex[0]] == 2):
                    for otherIndex in triangleWithoutIndex:
                        for edgeIndex in self.mesh.trianglesWithEdges[triangleindex]:
                            if self.mesh.edges[edgeIndex][0] in [index,otherIndex] and self.mesh.edges[edgeIndex][1] in [otherIndex,index]:
                                PotentialValue = self.fieldValues[index] + self.cost(self.density[triangleindex])*self.mesh.edgeLength[edgeIndex]
                                break
                        if PotentialValue < self.fieldValues[otherIndex]:
                            self.fieldValues[otherIndex] = PotentialValue
                            Visited = self.addInOrderByFieldValue(Visited,otherIndex)


        #Main loop
        LastValidated = Validated[0]
        NumVertices = len(StateMap)

        while(len(Validated) < NumVertices):

            #Validation of a vertex
            if(len(Visited) == 0):
                if(self.opt['debugging']):
                    print("Completion : ", len(Validated), "/", NumVertices)
                raise ValueError("The fast marching ended too fast... Is the Mesh connected ?")
            if(self.fieldValues[Visited[0]] == float('inf')):
                if(self.opt['debugging']):
                    print("Infinity value detected at vertex (",self.mesh.vertices[Visited[0]][0],",",self.mesh.vertices[Visited[0]][1], ")" )
                raise ValueError("Gradient computation impossible.")

            Validated.append(Visited[0])
            LastValidated = Visited[0]
            Visited = Visited[1:]
            StateMap[LastValidated] = 1

            #recomputation of the narrow band
            for [triangleindex, triangleWithoutIndex] in self.mesh.trianglesPerVertex[LastValidated]:

                if(StateMap[triangleWithoutIndex[1]] == 0):
                    StateMap[triangleWithoutIndex[1]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 0):
                    StateMap[triangleWithoutIndex[0]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 1 and StateMap[triangleWithoutIndex[1]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradConstrained(self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[LastValidated],
                                                                            self.fieldValues[triangleWithoutIndex[0]],
                                                                            self.fieldValues[LastValidated],
                                                                            self.cost(self.density[triangleindex]))
                    if PotentialValue < self.fieldValues[triangleWithoutIndex[1]]:
                        self.fieldValues[triangleWithoutIndex[1]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited, triangleWithoutIndex[1])


                elif(StateMap[triangleWithoutIndex[1]] == 1 and StateMap[triangleWithoutIndex[0]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradConstrained(self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[LastValidated],
                                                                            self.fieldValues[triangleWithoutIndex[1]],
                                                                            self.fieldValues[LastValidated],
                                                                            self.cost(self.density[triangleindex]))
                    if PotentialValue < self.fieldValues[triangleWithoutIndex[0]]:
                        self.fieldValues[triangleWithoutIndex[0]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited, triangleWithoutIndex[0])

                elif(StateMap[triangleWithoutIndex[1]] == 2 and StateMap[triangleWithoutIndex[0]] == 2):
                    for otherIndex in triangleWithoutIndex:
                        for edgeIndex in self.mesh.trianglesWithEdges[triangleindex]:
                            if self.mesh.edges[edgeIndex][0] in [LastValidated,otherIndex] and self.mesh.edges[edgeIndex][1] in [otherIndex,LastValidated]:
                                PotentialValue = self.fieldValues[LastValidated] + self.cost(self.density[triangleindex])*self.mesh.edgeLength[edgeIndex]
                                break
                        if PotentialValue < self.fieldValues[otherIndex]:
                            self.fieldValues[otherIndex] = PotentialValue
                            Visited = self.addInOrderByFieldValue(Visited,otherIndex)

    def computeFieldUnconstrainedDep1(self) -> None:
        """
        Computes the approximated solution of the eikonal equation and stores the approximation in 'fieldValues'.
        The method used is the FMT algorithm with constraints on the gradient and vertices at distance one or two are considered neighbour.
        i.e. constrained = False and NarrowBandDepth = 1 see 'opt'
        """
        StateMap = [0 for i in range(len(self.mesh.vertices))]

        Visited = []
        Validated = []
        self.fieldValues.setInfinity()

        #Dirichlet free exit condition
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            StateMap[index] = 1

        #Computation of the narrow band
        for index in Validated:
            for [triangleindex, triangleWithoutIndex] in self.mesh.trianglesPerVertex[index]:

                if(StateMap[triangleWithoutIndex[1]] == 0 and StateMap[triangleWithoutIndex[0]] == 1):
                    StateMap[triangleWithoutIndex[1]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 0 and StateMap[triangleWithoutIndex[1]] == 1):
                    StateMap[triangleWithoutIndex[0]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 1 and StateMap[triangleWithoutIndex[1]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradUnconstrained(self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[index],
                                                                            self.fieldValues[triangleWithoutIndex[0]],
                                                                            self.fieldValues[index],
                                                                            self.cost(self.density[triangleindex]))
                    if PotentialValue < self.fieldValues[triangleWithoutIndex[1]]:
                        self.fieldValues[triangleWithoutIndex[1]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[1])


                elif(StateMap[triangleWithoutIndex[1]] == 1 and StateMap[triangleWithoutIndex[0]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradUnconstrained(self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[index],
                                                                            self.fieldValues[triangleWithoutIndex[1]],
                                                                            self.fieldValues[index],
                                                                            self.cost(self.density[triangleindex]))
                    if PotentialValue < self.fieldValues[triangleWithoutIndex[0]]:
                        self.fieldValues[triangleWithoutIndex[0]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[0])

        #Main loop
        LastValidated = Validated[0]
        NumVertices = len(StateMap)

        while(len(Validated) < NumVertices):

            #Validation of a vertex
            if(len(Visited) == 0):
                if(self.opt['debugging']):
                    print("Completion : ", len(Validated), "/", NumVertices)
                raise ValueError("The fast marching ended too fast... Is the Mesh connected ?")
            if(self.fieldValues[Visited[0]] == float('inf')):
                if(self.opt['debugging']):
                    print("Infinity value detected at vertex (",self.mesh.vertices[Visited[0]][0],",",self.mesh.vertices[Visited[0]][1], ")" )
                raise ValueError("Gradient computation impossible.")

            Validated.append(Visited[0])
            LastValidated = Visited[0]
            Visited = Visited[1:]
            StateMap[LastValidated] = 1

            #Recomputation of the narrow band
            for [triangleindex, triangleWithoutIndex] in self.mesh.trianglesPerVertex[LastValidated]:

                if(StateMap[triangleWithoutIndex[1]] == 0 and StateMap[triangleWithoutIndex[0]] == 1):
                    StateMap[triangleWithoutIndex[1]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 0 and StateMap[triangleWithoutIndex[1]] == 1):
                    StateMap[triangleWithoutIndex[0]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 1 and StateMap[triangleWithoutIndex[1]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradUnconstrained(self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[LastValidated],
                                                                            self.fieldValues[triangleWithoutIndex[0]],
                                                                            self.fieldValues[LastValidated],
                                                                            self.cost(self.density[triangleindex]))
                    if PotentialValue < self.fieldValues[triangleWithoutIndex[1]]:
                        self.fieldValues[triangleWithoutIndex[1]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited, triangleWithoutIndex[1])



                elif(StateMap[triangleWithoutIndex[1]] == 1 and StateMap[triangleWithoutIndex[0]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradUnconstrained(self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[LastValidated],
                                                                            self.fieldValues[triangleWithoutIndex[1]],
                                                                            self.fieldValues[LastValidated],
                                                                            self.cost(self.density[triangleindex]))
                    if PotentialValue < self.fieldValues[triangleWithoutIndex[0]]:
                        self.fieldValues[triangleWithoutIndex[0]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited, triangleWithoutIndex[0])


    def computeFieldConstrainedDep1(self):
        """
        Computes the approximated solution of the eikonal equation and stores the approximation in 'fieldValues'.
        The method used is the FMT algorithm with constraints on the gradient and vertices at distance one or two are considered neighbour.
        i.e. constrained = True and NarrowBandDepth = 1 see 'opt'
        """
        StateMap = [0 for i in range(len(self.mesh.vertices))]

        Visited = []
        Validated = []
        self.fieldValues.setInfinity()

        #Dirichlet free exit condition
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            StateMap[index] = 1

        #Computation of the narrow band
        for index in Validated:
            for [triangleindex, triangleWithoutIndex] in self.mesh.trianglesPerVertex[index]:

                if(StateMap[triangleWithoutIndex[1]] == 0):
                    StateMap[triangleWithoutIndex[1]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 0):
                    StateMap[triangleWithoutIndex[0]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 1 and StateMap[triangleWithoutIndex[1]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradConstrained(self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[index],
                                                                            self.fieldValues[triangleWithoutIndex[0]],
                                                                            self.fieldValues[index],
                                                                            self.cost(self.density[triangleindex]))
                    if PotentialValue < self.fieldValues[triangleWithoutIndex[1]]:
                        self.fieldValues[triangleWithoutIndex[1]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[1])


                elif(StateMap[triangleWithoutIndex[1]] == 1 and StateMap[triangleWithoutIndex[0]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradConstrained(self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[index],
                                                                            self.fieldValues[triangleWithoutIndex[1]],
                                                                            self.fieldValues[index],
                                                                            self.cost(self.density[triangleindex]))
                    if PotentialValue < self.fieldValues[triangleWithoutIndex[0]]:
                        self.fieldValues[triangleWithoutIndex[0]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited,triangleWithoutIndex[0])

                elif(StateMap[triangleWithoutIndex[1]] == 2 and StateMap[triangleWithoutIndex[0]] == 2):
                    PotentialValue = self.fieldValues[index] + self.cost(self.density[triangleindex])*EikoSolver.computeHeightLength(self.mesh.vertices[index],
                                                                                                    self.mesh.vertices[triangleWithoutIndex[0]],
                                                                                                    self.mesh.vertices[triangleWithoutIndex[1]])

                    for i in range(2):
                        if PotentialValue < self.fieldValues[triangleWithoutIndex[i]]:
                            self.fieldValues[triangleWithoutIndex[i]] = PotentialValue
                            Visited = self.addInOrderByFieldValue(Visited, triangleWithoutIndex[i])

        #Main loop
        LastValidated = Validated[0]
        NumVertices = len(StateMap)

        while(len(Validated) < NumVertices):

            #Validation of a vertex
            if(len(Visited) == 0):
                if(self.opt['debugging']):
                    print("Completion : ", len(Validated), "/", NumVertices)
                raise ValueError("The fast marching ended too fast... Is the Mesh connected ?")
            if(self.fieldValues[Visited[0]] == float('inf')):
                if(self.opt['debugging']):
                    print("Infinity value detected at vertex (",self.mesh.vertices[Visited[0]][0],",",self.mesh.vertices[Visited[0]][1], ")" )
                raise ValueError("Gradient computation impossible.")

            Validated.append(Visited[0])
            LastValidated = Visited[0]
            Visited = Visited[1:]
            StateMap[LastValidated] = 1

            #Recomputation of the narrow band
            for [triangleindex, triangleWithoutIndex] in self.mesh.trianglesPerVertex[LastValidated]:

                if(StateMap[triangleWithoutIndex[1]] == 0):
                    StateMap[triangleWithoutIndex[1]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 0):
                    StateMap[triangleWithoutIndex[0]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 1 and StateMap[triangleWithoutIndex[1]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradConstrained(self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[LastValidated],
                                                                            self.fieldValues[triangleWithoutIndex[0]],
                                                                            self.fieldValues[LastValidated],
                                                                            self.cost(self.density[triangleindex]))
                    if PotentialValue < self.fieldValues[triangleWithoutIndex[1]]:
                        self.fieldValues[triangleWithoutIndex[1]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited, triangleWithoutIndex[1])


                elif(StateMap[triangleWithoutIndex[1]] == 1 and StateMap[triangleWithoutIndex[0]] == 2):
                    PotentialValue = EikoSolver.computeHeightFromGradConstrained(self.mesh.vertices[triangleWithoutIndex[0]],
                                                                            self.mesh.vertices[triangleWithoutIndex[1]],
                                                                            self.mesh.vertices[LastValidated],
                                                                            self.fieldValues[triangleWithoutIndex[1]],
                                                                            self.fieldValues[LastValidated],
                                                                            self.cost(self.density[triangleindex]))
                    if PotentialValue < self.fieldValues[triangleWithoutIndex[0]]:
                        self.fieldValues[triangleWithoutIndex[0]] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited, triangleWithoutIndex[0])

                elif(StateMap[triangleWithoutIndex[1]] == 2 and StateMap[triangleWithoutIndex[0]] == 2):
                    PotentialValue = self.fieldValues[LastValidated] + self.cost(self.density[triangleindex])*EikoSolver.computeHeightLength(self.mesh.vertices[LastValidated],
                                                                                                    self.mesh.vertices[triangleWithoutIndex[0]],
                                                                                                    self.mesh.vertices[triangleWithoutIndex[1]])

                    for i in range(2):
                        if PotentialValue < self.fieldValues[triangleWithoutIndex[i]]:
                            self.fieldValues[triangleWithoutIndex[i]] = PotentialValue
                            Visited = self.addInOrderByFieldValue(Visited, triangleWithoutIndex[i])


    def computeFieldByEdges(self) -> None:
        """
        Computes a numerical approximation of the solution to the eikonal equation using a FME algorithm.
        """

        StateMap = [0 for i in range(len(self.mesh.vertices))]

        Visited = []
        Validated = []
        self.fieldValues.setInfinity()

        #Dirichlet free exit condition
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            StateMap[index] = 1

        #Computation of the narrow band
        for index in Validated:
            for [triangleindex, triangleWithoutIndex] in self.mesh.trianglesPerVertex[index]:
                for otherIndex in triangleWithoutIndex:
                    if(StateMap[otherIndex] == 0):
                        StateMap[otherIndex] = 2
                    for edgeIndex in self.mesh.trianglesWithEdges[triangleindex]:
                        if self.mesh.edges[edgeIndex][0] in [index,otherIndex] and self.mesh.edges[edgeIndex][1] in [otherIndex,index]:
                            PotentialValue = self.fieldValues[index] + self.cost(self.density[triangleindex])*self.mesh.edgeLength[edgeIndex]
                            break
                    if PotentialValue < self.fieldValues[otherIndex]:
                        self.fieldValues[otherIndex] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited,otherIndex)

        #Main loop
        LastValidated = Validated[0]
        NumVertices = len(StateMap)

        while(len(Validated) < NumVertices):

            #Validation of a vertex
            if(len(Visited) == 0):
                if(self.opt['debugging']):
                    print("Completion : ", len(Validated), "/", NumVertices)
                raise ValueError("The fast marching ended too fast... Is the Mesh connected ?")
            if(self.fieldValues[Visited[0]] == float('inf')):
                if(self.opt['debugging']):
                    print("Infinity value detected at vertex (",self.mesh.vertices[Visited[0]][0],",",self.mesh.vertices[Visited[0]][1], ")" )
                raise ValueError("Gradient computation impossible.")

            Validated.append(Visited[0])
            LastValidated = Visited[0]
            Visited = Visited[1:]
            StateMap[LastValidated] = 1

            #Recomputation of the narrow band
            for [triangleindex, triangleWithoutIndex] in self.mesh.trianglesPerVertex[LastValidated]:

                for otherIndex in triangleWithoutIndex:
                    if(StateMap[otherIndex] == 0):
                        StateMap[otherIndex] = 2
                    for edgeIndex in self.mesh.trianglesWithEdges[triangleindex]:
                        if self.mesh.edges[edgeIndex][0] in [LastValidated,otherIndex] and self.mesh.edges[edgeIndex][1] in [otherIndex,LastValidated]:
                            PotentialValue = self.fieldValues[LastValidated] + self.cost(self.density[triangleindex])*self.mesh.edgeLength[edgeIndex]
                            break
                    if PotentialValue < self.fieldValues[otherIndex]:
                        self.fieldValues[otherIndex] = PotentialValue
                        Visited = self.addInOrderByFieldValue(Visited,otherIndex)

    def computeField(self) -> None:
        """
        Computes the approximated solution of the eikonal equation and stores the approximation in 'fieldValues'.
        The method used depends on the options set in the 'opt' dictionary.
        """
        if(self.opt['method'] == "FMT"):
            if(self.opt['constrained'] and self.opt['NarrowBandDepth'] == 2):
                self.computeFieldConstrainedDep2()
            elif( not self.opt['constrained'] and self.opt['NarrowBandDepth'] == 2):
                self.computeFieldUnconstrainedDep2()
            elif( self.opt['constrained'] and self.opt['NarrowBandDepth'] == 1):
                self.computeFieldConstrainedDep1()
            elif( not self.opt['constrained'] and self.opt['NarrowBandDepth'] == 1):
                self.computeFieldUnconstrainedDep1()
        elif(self.opt['method'] == "FME"):
            self.computeFieldByEdges()

    def showNarrowBandAfterStep(self, n:int) -> None:
        """
        Computes the first `n` vertices of the approximation, then, displays the narrow band.

        Args:
            n (int): number of vertices to compute.
        """
        NotVisited = [i for i in range(len(self.mesh.vertices))]

        Visited = []
        Validated = []
        self.fieldValues.setInfinity()

        #Dirichlet free exit condition
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            NotVisited.remove(index)

        #Computation of the narrow band
        for index in Validated:
            for triangleindex, triangle in enumerate(self.mesh.triangles):
                if index in triangle:
                    count = 0
                    selected = -1
                    ValidatedPointindex = -1
                    for j,i in enumerate(triangle):
                        if(i in Validated):
                            count += 1
                            offsetValidatedPoint = j
                        elif i in NotVisited:
                            selected = i
                            offset = j
                        elif i in Visited:
                            selected = i
                            offset = j
                    if(count == 2 and selected != -1):
                        if(selected not in Visited):
                            Visited.append(selected)
                            NotVisited.remove(selected)

                        if(self.opt['constrained']):
                            PotentialValue = EikoSolver.computeHeightFromGrad(self.mesh.vertices[triangle[offset]],
                                                                                self.mesh.vertices[triangle[(offset + 1)%3]],
                                                                                self.mesh.vertices[triangle[(offset + 2)%3]],
                                                                                self.fieldValues[triangle[(offset + 1)%3]],
                                                                                self.fieldValues[triangle[(offset + 2)%3]],
                                                                                self.cost(self.density[triangleindex]))
                        else:
                            PotentialValue = EikoSolver.computeHeightFromGrad2(self.mesh.vertices[triangle[offset]],
                                                                                self.mesh.vertices[triangle[(offset + 1)%3]],
                                                                                self.mesh.vertices[triangle[(offset + 2)%3]],
                                                                                self.fieldValues[triangle[(offset + 1)%3]],
                                                                                self.fieldValues[triangle[(offset + 2)%3]],
                                                                                self.cost(self.density[triangleindex]))
                        if PotentialValue < self.fieldValues[selected]:
                            self.fieldValues[selected] = PotentialValue



                    elif(self.opt['NarrowBandDepth'] == 2 and count == 1):

                        for i in triangle:
                            if(i in NotVisited):
                                Visited.append(i)
                                NotVisited.remove(i)

                            if(i not in Validated):
                                PotentialValue = self.fieldValues[triangle[offsetValidatedPoint]] + self.cost(self.density[triangleindex])*EikoSolver.computeHeightLength(self.mesh.vertices[triangle[offsetValidatedPoint]],
                                                                                                                self.mesh.vertices[triangle[(offsetValidatedPoint + 1)%3]],
                                                                                                                self.mesh.vertices[triangle[(offsetValidatedPoint + 2)%3]])
                                if PotentialValue < self.fieldValues[i]:
                                    self.fieldValues[i] = PotentialValue

        #Main loop
        LastValidated = Validated[0]
        for vide in range(n):

            #Recomputation of the narrow band
            for triangleindex,triangle in enumerate(self.mesh.triangles):
                if(LastValidated in triangle):
                    count = 0
                    selected = -1
                    for j,i in enumerate(triangle):
                        if(i in Validated):
                            count += 1
                            offsetValidatedPoint = j
                        elif i in NotVisited:
                            selected = i
                            offset = j
                        elif i in Visited:
                            selected = i
                            offset = j
                    if(count == 2 and selected != -1):
                        if(selected not in Visited):
                            Visited.append(selected)
                            NotVisited.remove(selected)

                        if(self.opt['constrained']):
                            PotentialValue = EikoSolver.computeHeightFromGrad(self.mesh.vertices[triangle[offset]],
                                                                                self.mesh.vertices[triangle[(offset + 1)%3]],
                                                                                self.mesh.vertices[triangle[(offset + 2)%3]],
                                                                                self.fieldValues[triangle[(offset + 1)%3]],
                                                                                self.fieldValues[triangle[(offset + 2)%3]],
                                                                                self.cost(self.density[triangleindex]))
                        else:
                            PotentialValue = EikoSolver.computeHeightFromGrad2(self.mesh.vertices[triangle[offset]],
                                                                                self.mesh.vertices[triangle[(offset + 1)%3]],
                                                                                self.mesh.vertices[triangle[(offset + 2)%3]],
                                                                                self.fieldValues[triangle[(offset + 1)%3]],
                                                                                self.fieldValues[triangle[(offset + 2)%3]],
                                                                                self.cost(self.density[triangleindex]))
                        if PotentialValue < self.fieldValues[selected]:
                            self.fieldValues[selected] = PotentialValue



                    elif(self.opt['NarrowBandDepth'] == 2 and count == 1):

                        for i in triangle:
                            if(i in NotVisited):
                                Visited.append(i)
                                NotVisited.remove(i)

                            if(i not in Validated):
                                PotentialValue = self.fieldValues[triangle[offsetValidatedPoint]] + self.cost(self.density[triangleindex])*EikoSolver.computeHeightLength(self.mesh.vertices[triangle[offsetValidatedPoint]],
                                                                                                                self.mesh.vertices[triangle[(offsetValidatedPoint + 1)%3]],
                                                                                                                self.mesh.vertices[triangle[(offsetValidatedPoint + 2)%3]])
                                if PotentialValue < self.fieldValues[i]:
                                    self.fieldValues[i] = PotentialValue

            #Validation:
            minIndex = Visited[0]
            minValue = self.fieldValues[minIndex]
            for index in Visited:
                if minValue > self.fieldValues[index]:
                    minIndex = index
                    minValue = self.fieldValues[index]

            if(minValue == float('inf')):
                if(self.opt['debugging']):
                    print("Infinity value detected at vertex (",self.mesh.vertices[Visited[0]][0],",",self.mesh.vertices[Visited[0]][1], ")" )
                raise ValueError("Gradient computation impossible.")

            Validated.append(minIndex)
            LastValidated = minIndex
            Visited.remove(minIndex)

        Visited.append(minIndex)
        Validated.remove(minIndex)

        for index in Validated:
            self.fieldValues[index] = 0.5
        for index in Visited:
            self.fieldValues[index] = 0
        for index in NotVisited:
            self.fieldValues[index] = 1

        self.fieldValues.show(grid=True)

    def addInOrderByFieldValue(self,L:List[int],index:int) -> List[int]:
        """
        Adds the index given as an argument to the given list. Returns the new list. The index is inserted such that the `fieldValues` corresponding to the vertices of the indices are increasing.

        Example:
            If we have the following lists::

                MyMesh = Mesh()
                MyMesh.vertices = [[0,0],[0,1],[1,1],[1,0]]
                MySolver = EikoSolver(MyMesh)
                MySolver.fieldValues.values = [0.5,4.5,2.0,3.0]
                MyList = [0,3]

            Then::

                >>>MySolver.addInOrderByFieldValue(MyList,1)
                [0,3,1]
                >>>MySolver.addInOrderByFieldValue(MyList,2)
                [0,2,3]

        Note:
            If an index is already present in the list, its position in the list is recomputed.

        Args:
            L (List[int]): the list in which the new index will be inserted.
            index (int): the index to insert.

        Returns:
            List[int]: the new list containing the index.
        """
        try:
            L.remove(index)
        except ValueError:
            pass # or scream: thing not in some_list!

        if(len(L) == 0):
            return [index]
        rank = 0
        while(self.fieldValues[index] > self.fieldValues[L[rank]]):
            rank += 1
            if rank >= len(L):
                break

        return L[:rank] + [index] + L[rank:]


    @staticmethod
    def computeHeightFromGradUnconstrained(C:List[float],B:List[float],A:List[float],Vb:float,Va:float,P:float) -> float:
        """
        Computes and returns the value of Vc such that, if we denote by :math:`\\Phi_{ABC}(Va,Vb,Vc)` the unique affine function :math:`F` of :math:`\\mathbb{R}^2` such that
        .. math::

            F(A) = Va, \\; F(B) = Vb, \\; F(C) = Vc

        then we have
        .. math::

            \\left| \\nabla\\Phi_{ABC}(Va,Vb,Vc) \\right| = P.

        Args:
            C (List[float]): a point of the triangle.
            B (List[float]): a point of the triangle.
            A (List[float]): a point of the triangle.
            Vb (float): the value of F at B.
            Va (float): the value of F at A.
            P (float): the norm of the gradient.

        Returns:
            float: the value Vc.
        """
        if(Vb == float('inf') or Va == float('inf') or P < 0):
            if self.opt["debugging"]:
                print("Vb = ", Vb, ", Va = ", Va, ", P = ",P)
            raise ValueError("Invalid values for the gradient constraint")
        AB2 = (A[0] - B[0])*(A[0] - B[0]) + (A[1] - B[1])*(A[1] - B[1])


        if(AB2*(P**2) - ( (Va - Vb)**2) <0):
            if self.opt["debugging"]:
                print("Vb = ", Vb, ", Va = ", Va, ", P = ",P)
            print("Warning - Given values are incompatible with the gradient constraint")
            AC2 = (A[0] - C[0])*(A[0] - C[0]) + (A[1] - C[1])*(A[1] - C[1])
            return Va + np.sqrt(AC2)*P


        detCACB = abs((C[0] - A[0])*(C[1] - B[1]) - (C[1] - A[1])*(C[0] - B[0]))
        CBxCA = (B[0] - C[0])*(A[0] - C[0]) + (B[1] - C[1])*(A[1] - C[1])

        if(Va > Vb):
            BC2 = (C[0] - B[0])*(C[0] - B[0]) + (C[1] - B[1])*(C[1] - B[1])
            return Vb + (Va-Vb)*(BC2 - CBxCA)/AB2 + detCACB*np.sqrt(AB2*(P**2) - ( (Va - Vb)**2) )/AB2
        else:
            AC2 = (A[0] - C[0])*(A[0] - C[0]) + (A[1] - C[1])*(A[1] - C[1])
            return Va + (Vb-Va)*(AC2 - CBxCA)/AB2 + detCACB*np.sqrt(AB2*(P**2) - ( (Vb - Va)**2) )/AB2

    def computeHeightFromGradConstrained(A:List[float],B:List[float],C:List[float],Vb:float,Vc:float,P:float) -> float:
        """
        Computes and returns the value of Va such that, if we denote by :math:`\\Phi_{ABC}(Va,Vb,Vc)` the unique affine function :math:`F` of :math:`\\mathbb{R}^2` such that
        .. math::

            F(A) = Va, \\; F(B) = Vb, \\; F(C) = Vc

        then we have
        .. math::

            \\left| \\nabla\\Phi_{ABC}(Va,Vb,Vc) \\right| = P.

        If the gradient found is not inside the triangle, Va is taken as the smallest value found by following an edge of the triangle.

        Args:
            A (List[float]): a point of the triangle.
            B (List[float]): a point of the triangle.
            C (List[float]): a point of the triangle.
            Vb (float): the value of F at B.
            Vc (float): the value of F at C.
            P (float): the norm of the gradient.

        Returns:
            float: the value Va.
        """
        if(Vb == float('inf') or Vc == float('inf') or P < 0):
            if self.opt["debugging"]:
                print("Vb = ", Vb, ", Vc = ", Vc, ", P = ",P)
            raise ValueError("Invalid values for the gradient constraint")

        if( Vb > Vc):
            s1 = Vc
            s2 = list(C)
            C = B
            Vc = Vb
            B = s2
            Vb = s1

        AB2 = (A[0] - B[0])*(A[0] - B[0]) + (A[1] - B[1])*(A[1] - B[1])
        ABxBC = (B[0] - C[0])*(A[0] - B[0]) + (B[1] - C[1])*(A[1] - B[1])

        if(P*ABxBC > (Vc-Vb)*np.sqrt(AB2)):
            return Vb + np.sqrt(AB2)*P

        AC2 = (A[0] - C[0])*(A[0] - C[0]) + (A[1] - C[1])*(A[1] - C[1])
        CBxCA = (B[0] - C[0])*(A[0] - C[0]) + (B[1] - C[1])*(A[1] - C[1])

        if(P*CBxCA < (Vc-Vb)*np.sqrt(AC2)):
            return Vc + np.sqrt(AC2)*P

        BC2 = (B[0] - C[0])*(B[0] - C[0]) + (B[1] - C[1])*(B[1] - C[1])

        if(BC2*P**2 - ( (Vc - Vb)**2) <0):
            if self.opt["debugging"]:
                print("Vb = ", Vb, ", Vc = ", Vc, ", P = ",P)
            print("Warning - Given values are incompatible with the gradient constraint")
            return Vc + np.sqrt(AC2)*P

        detABCB = abs((B[0] - A[0])*(B[1] - C[1]) - (B[1] - A[1])*(B[0] - C[0]))

        return Vb - (Vc-Vb)*(ABxBC)/BC2 + detABCB*np.sqrt(BC2*P**2 - ( (Vb - Vc)**2) )/BC2

    def computeHeightLength(B:float,C:float,A:float) -> float:
        """
        Computes and returns the length of the height of the triangle ABC passing through B.

        Args:
            B (List[float]): a point of the triangle.
            C (List[float]): a point of the triangle.
            A (List[float]): a point of the triangle.

        Returns:
            float: the length of the height.
        """
        detABAC = abs((A[0] - C[0])*(A[1] - B[1]) - (A[1] - C[1])*(A[0] - B[0]))
        AC2 = (A[0] - C[0])*(A[0] - C[0]) + (A[1] - C[1])*(A[1] - C[1])
        return detABAC/np.sqrt(AC2)

    def dist(x0:float,y0:float,x1:float,y1:float) -> float:
        """
        Computes and returns the euclidian distance between (x0,y0) and (x1,y1).
        """
        return(np.sqrt((x1-x0)**2 + (y1-y0)**2))
