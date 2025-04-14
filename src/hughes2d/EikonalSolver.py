#from Mesh2D import *
from hughes2d.Mesh2D import *
#import plotly
import random as alea


class EikoSolver(object):

    def __init__(self, Mesh, DensityMap = [], costFunction = (lambda x: 1+x), opt=dict()):
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
            self.opt['constrained'] = False

        if("debugging" not in self.opt.keys()):
            self.opt['debugging'] = False

        self.cost = costFunction

    def updateDensity(self, density):
        self.density = density

    def computeFieldUnconstrainedDep2(self):
        StateMap = [0 for i in range(len(self.mesh.vertices))]
        #print(len(self.mesh.vertices))
        VisitedDict = dict()
        Visited = []
        TriangleOfVisit = []
        Validated = []
        self.fieldValues.setInfinity()


        #On ajoute les dirichlet sortie libre
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            StateMap[index] = 1


        #Selection de la liste de la nouvelle generation a traiter

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


        #Boucle principale
        LastValidated = Validated[0]
        NumVertices = len(StateMap)



        while(len(Validated) < NumVertices):

            #Valider le plus petit:
            if(len(Visited) == 0):
                print("Fin prematuree.... Num Validated : ", len(Validated))
            if(Visited[0] == float('inf')):
                print("Impossible de calculer le gradient")
                print(Visited)
                print(NotVisited)
                self.fieldValues[minIndex] = 1000000000

            Validated.append(Visited[0])
            LastValidated = Visited[0]
            Visited = Visited[1:]
            StateMap[LastValidated] = 1

            """
            print("Value validated : ", self.fieldValues[LastValidated])
            print("True Value : ", 10 - self.mesh.vertices[LastValidated][0])
            print("Point : ", self.mesh.vertices[LastValidated])
            """

            #Ajout des nouveaux consider
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

    def computeFieldConstrainedDep2(self):
        StateMap = [0 for i in range(len(self.mesh.vertices))]
        #print(len(self.mesh.vertices))
        Visited = []
        Validated = []
        self.fieldValues.setInfinity()


        #On ajoute les dirichlet sortie libre
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            StateMap[index] = 1


        #Selection de la liste de la nouvelle generation a traiter

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



        #print("Premier element : ", Validated[0])


        #Boucle principale
        LastValidated = Validated[0]
        NumVertices = len(StateMap)



        while(len(Validated) < NumVertices):

            #Valider le plus petit:
            if(len(Visited) == 0):
                print("Fin prematuree.... Num Validated : ", len(Validated))
            if(Visited[0] == float('inf')):
                print("Impossible de calculer le gradient")
                print(Visited)
                print(NotVisited)
                self.fieldValues[minIndex] = 1000000000

            Validated.append(Visited[0])
            LastValidated = Visited[0]
            Visited = Visited[1:]
            StateMap[LastValidated] = 1


            #Ajout des nouveaux consider
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

    def computeFieldUnconstrainedDep1(self):
        StateMap = [0 for i in range(len(self.mesh.vertices))]
        #print(len(self.mesh.vertices))
        Visited = []
        Validated = []
        self.fieldValues.setInfinity()

        #On ajoute les dirichlet sortie libre
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            StateMap[index] = 1


        #Selection de la liste de la nouvelle generation a traiter

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
                        if(self.fieldValues[Visited[0]] > self.fieldValues[Visited[1]]):
                            raise ValueError("ICI !!!!")

        #Boucle principale
        LastValidated = Validated[0]
        NumVertices = len(StateMap)



        while(len(Validated) < NumVertices):

            #Valider le plus petit:
            if(len(Visited) == 0):
                print("Fin prematuree.... Num Validated : ", len(Validated))
            if(Visited[0] == float('inf')):
                print("Impossible de calculer le gradient")
                print(Visited)
                print(NotVisited)
                self.fieldValues[minIndex] = 1000000000

            Validated.append(Visited[0])
            LastValidated = Visited[0]
            Visited = Visited[1:]
            StateMap[LastValidated] = 1


            #Ajout des nouveaux consider
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
        StateMap = [0 for i in range(len(self.mesh.vertices))]
        #print(len(self.mesh.vertices))
        Visited = []
        Validated = []
        self.fieldValues.setInfinity()


        #On ajoute les dirichlet sortie libre
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            StateMap[index] = 1


        #Selection de la liste de la nouvelle generation a traiter

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



        #print("Premier element : ", Validated[0])


        #Boucle principale
        LastValidated = Validated[0]
        NumVertices = len(StateMap)



        while(len(Validated) < NumVertices):

            #Valider le plus petit:
            if(len(Visited) == 0):
                print("Fin prematuree.... Num Validated : ", len(Validated))
            if(Visited[0] == float('inf')):
                print("Impossible de calculer le gradient")
                print(Visited)
                print(NotVisited)
                self.fieldValues[minIndex] = 1000000000

            Validated.append(Visited[0])
            LastValidated = Visited[0]
            Visited = Visited[1:]
            StateMap[LastValidated] = 1


            #Ajout des nouveaux consider
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
        """StateMap = [0 for i in range(len(self.mesh.vertices))]
        #print(len(self.mesh.vertices))
        Visited = []
        Validated = []
        self.fieldValues.setInfinity()

        #On ajoute les dirichlet sortie libre
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            StateMap[index] = 1


        #Selection de la liste de la nouvelle generation a traiter

        for index in Validated:
            for [triangleindex, triangleWithoutIndex] in self.mesh.trianglesPerVertex[index]:

                if(StateMap[triangleWithoutIndex[1]] == 0 and StateMap[triangleWithoutIndex[0]] == 1):
                    StateMap[triangleWithoutIndex[1]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 0 and StateMap[triangleWithoutIndex[1]] == 1):
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


        #Boucle principale
        LastValidated = Validated[0]
        NumVertices = len(StateMap)



        while(len(Validated) < NumVertices):

            #Valider le plus petit:
            if(len(Visited) == 0):
                print("Fin prematuree.... Num Validated : ", len(Validated))
            if(Visited[0] == float('inf')):
                print("Impossible de calculer le gradient")
                print(Visited)
                print(NotVisited)
                self.fieldValues[minIndex] = 1000000000

            Validated.append(Visited[0])
            LastValidated = Visited[0]
            Visited = Visited[1:]
            StateMap[LastValidated] = 1


            #Ajout des nouveaux consider
            for [triangleindex, triangleWithoutIndex] in self.mesh.trianglesPerVertex[LastValidated]:

                if(StateMap[triangleWithoutIndex[1]] == 0 and StateMap[triangleWithoutIndex[0]] == 1):
                    StateMap[triangleWithoutIndex[1]] = 2

                if(StateMap[triangleWithoutIndex[0]] == 0 and StateMap[triangleWithoutIndex[1]] == 1):
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
                        Visited = self.addInOrderByFieldValue(Visited, triangleWithoutIndex[0])"""

    def computeFieldByEdges(self):
        """
        Computes a numerical approximation of the solution to the eikonal equation using a FME algorithm.
        """

        StateMap = [0 for i in range(len(self.mesh.vertices))]
        #print(len(self.mesh.vertices))
        Visited = []
        Validated = []
        self.fieldValues.setInfinity()


        #On ajoute les dirichlet sortie libre
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            StateMap[index] = 1


        #Selection de la liste de la nouvelle generation a traiter

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

        #Boucle principale
        LastValidated = Validated[0]
        NumVertices = len(StateMap)


        while(len(Validated) < NumVertices):

            #Valider le plus petit:
            if(len(Visited) == 0):
                print("Fin prematuree.... Num Validated : ", len(Validated))
            if(Visited[0] == float('inf')):
                print("Impossible de calculer le gradient")
                print(Visited)
                print(NotVisited)
                self.fieldValues[minIndex] = 1000000000

            Validated.append(Visited[0])
            LastValidated = Visited[0]
            Visited = Visited[1:]
            StateMap[LastValidated] = 1


            #Ajout des nouveaux consider
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



    def computeField(self):
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

    def showNarrowBandAfterStep(self, n):
        NotVisited = [i for i in range(len(self.mesh.vertices))]
        #print(len(self.mesh.vertices))
        Visited = []
        Validated = []
        self.fieldValues.setInfinity()

        #On ajoute les dirichlet sortie libre
        for index in self.mesh.exitVertices:
            self.fieldValues[index] = 0
            Validated.append(index)
            NotVisited.remove(index)


        #Selection de la liste de la nouvelle generation a traiter

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
                                """PotentialValue = self.fieldValues[ValidatedPointindex] + self.cost(self.density[triangleindex])*EikoSolver.dist(self.mesh.vertices[triangle[offsetValidatedPoint]][0],
                                                                                                                self.mesh.vertices[triangle[offsetValidatedPoint]][1],
                                                                                                                self.mesh.vertices[triangle[i]][0],
                                                                                                                self.mesh.vertices[triangle[i]][1])"""

                                PotentialValue = self.fieldValues[triangle[offsetValidatedPoint]] + self.cost(self.density[triangleindex])*EikoSolver.computeHeightLength(self.mesh.vertices[triangle[offsetValidatedPoint]],
                                                                                                                self.mesh.vertices[triangle[(offsetValidatedPoint + 1)%3]],
                                                                                                                self.mesh.vertices[triangle[(offsetValidatedPoint + 2)%3]])
                                if PotentialValue < self.fieldValues[i]:
                                    self.fieldValues[i] = PotentialValue



        #print("Premier element : ", Validated[0])


        #Boucle principale
        LastValidated = Validated[0]
        for vide in range(n):

            #Ajout des nouveaux consider
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
                                """PotentialValue = self.fieldValues[ValidatedPointindex] + self.cost(self.density[triangleindex])*EikoSolver.dist(self.mesh.vertices[triangle[offsetValidatedPoint]][0],
                                                                                                                self.mesh.vertices[triangle[offsetValidatedPoint]][1],
                                                                                                                self.mesh.vertices[triangle[i]][0],
                                                                                                                self.mesh.vertices[triangle[i]][1])"""

                                PotentialValue = self.fieldValues[triangle[offsetValidatedPoint]] + self.cost(self.density[triangleindex])*EikoSolver.computeHeightLength(self.mesh.vertices[triangle[offsetValidatedPoint]],
                                                                                                                self.mesh.vertices[triangle[(offsetValidatedPoint + 1)%3]],
                                                                                                                self.mesh.vertices[triangle[(offsetValidatedPoint + 2)%3]])
                                if PotentialValue < self.fieldValues[i]:
                                    self.fieldValues[i] = PotentialValue

            #Valider le plus petit:
            minIndex = Visited[0]
            minValue = self.fieldValues[minIndex]
            for index in Visited:
                if minValue > self.fieldValues[index]:
                    minIndex = index
                    minValue = self.fieldValues[index]

            if(minValue == float('inf')):
                print("Impossible de calculer le gradient")
                print(Visited)
                print(NotVisited)
                self.fieldValues[minIndex] = 1000000000

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

    def addInOrderByFieldValue(self,L,index):
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
    def computeHeightFromGradUnconstrained(C,B,A,Vb,Va,P): #Calcul de u(A) pour que grad du triangle ait la norme prescrite par P
        if(Vb == float('inf') or Va == float('inf') or P < 0):
            print("ERREUR ICI !!!!")
        AB2 = (A[0] - B[0])*(A[0] - B[0]) + (A[1] - B[1])*(A[1] - B[1])


        if(AB2*(P**2) - ( (Va - Vb)**2) <0):
            print("Gros probleme, ca ne devrait pas être possible")
            AC2 = (A[0] - C[0])*(A[0] - C[0]) + (A[1] - C[1])*(A[1] - C[1])
            return Va + np.sqrt(AC2)*P #valeur en longeant AC ?
            #return float('inf')

        detCACB = abs((C[0] - A[0])*(C[1] - B[1]) - (C[1] - A[1])*(C[0] - B[0]))
        CBxCA = (B[0] - C[0])*(A[0] - C[0]) + (B[1] - C[1])*(A[1] - C[1])

        if(Va > Vb):
            BC2 = (C[0] - B[0])*(C[0] - B[0]) + (C[1] - B[1])*(C[1] - B[1])
            return Vb + (Va-Vb)*(BC2 - CBxCA)/AB2 + detCACB*np.sqrt(AB2*(P**2) - ( (Va - Vb)**2) )/AB2
        else:
            AC2 = (A[0] - C[0])*(A[0] - C[0]) + (A[1] - C[1])*(A[1] - C[1])
            return Va + (Vb-Va)*(AC2 - CBxCA)/AB2 + detCACB*np.sqrt(AB2*(P**2) - ( (Vb - Va)**2) )/AB2

    def computeHeightFromGradConstrained(A,B,C,Vb,Vc,P): #Calcul de u(A) pour que grad du triangle ait la norme prescrite par P
        if(Vb == float('inf') or Vc == float('inf') or P < 0):
            print("ERREUR ICI !!!!")

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
            print("ERREUR cas non traité !!!!!!!!!!!!")
            return Vc + np.sqrt(AC2)*P #valeur en longeant AC ?
            #return float('inf')

        detABCB = abs((B[0] - A[0])*(B[1] - C[1]) - (B[1] - A[1])*(B[0] - C[0]))

        return Vb - (Vc-Vb)*(ABxBC)/BC2 + detABCB*np.sqrt(BC2*P**2 - ( (Vb - Vc)**2) )/BC2

    def computeHeightLength(B,C,A): #calcul de la longueur de la hauteur issue de B
        detABAC = abs((A[0] - C[0])*(A[1] - B[1]) - (A[1] - C[1])*(A[0] - B[0]))
        AC2 = (A[0] - C[0])*(A[0] - C[0]) + (A[1] - C[1])*(A[1] - C[1])
        return detABAC/np.sqrt(AC2)

    def dist(x0,y0,x1,y1):
        return(np.sqrt((x1-x0)**2 + (y1-y0)**2))
