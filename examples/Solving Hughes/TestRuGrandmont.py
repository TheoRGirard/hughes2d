from hughes2d import *
import random as alea
#import plotly



Exits = [[[2,3.5],[1.5,3.5]],[[5.25,4.25],[5.75,4.25]],[[9.75,4.25],[9.75,4.75]]]

ListOuterWall2 = [[8.5,10.75],[9,10],[8.25,9.5],[9.25,8.25],[10,8.75],[10.25,8.5],
                    [9.75,7.75],[10.25,7.25],[10.75,7.75],[12.25,5.5],[13,4],[12,3.75],[10.25,3.50],
                    [8.25,3.25],[6.25,3.25],[6.25,4.25],[5.75,4.25],[5.25,4.25],[5.25,3.25],
                    [4.75,3.25],[2,3.5],[1.5,3.5],[1,3.5],[2,5.25],[4,8.25],[5.25,7.25],
                    [6.25,8.5],[6.25,9.5],[5.5,10],[6,10.75],[6.25,10.5],[6.5,10.25],
                    [7,10],[7.5,10],[8,10.25],[8.25,10.5]]

Entries = [[[6.25,10.5],[6.5,10.25]],[[8,10.25],[8.25,10.5]]]

InnerWall1 = [[7,9.5],[7.75,9.5],[8,9.25],[8,8.75],[7.75,8.5],[7,8.5]]

InnerWall2 = [[6.25,8],[8.75,8],[9.25,7.25],[9.5,7.5],[9.75,7.25],[9.25,6.75],[10.25,5.5],[10.25,5],[10,4.75],[9.75,4.75],[9.75,4.25],[8.25,4.25],
                [8.25,4.5],[5,4.5],[4.25,4.75],[4.25,5.25],[5.75,7.25]]

for i in range(len(ListOuterWall2)):
    ListOuterWall2[i][1] = - ListOuterWall2[i][1]

for i in range(len(Exits)):
    Exits[i][0][1] = - Exits[i][0][1]
    Exits[i][1][1] = - Exits[i][1][1]

for i in range(len(InnerWall1)):
    InnerWall1[i][1] = - InnerWall1[i][1]

for i in range(len(InnerWall2)):
    InnerWall2[i][1] = - InnerWall2[i][1]

#Construction du domaine--------------------------------
MyDomain = Mesh2D.NonConvexDomain(ListOuterWall2)
MyDomain.addWall(InnerWall1,cycle=True)
MyDomain.addWall(InnerWall2,cycle=True)
MyDomain.addExits(Exits)
MyDomain.show()

#construction du Mesh--------------------------------------
MyMesh = Mesh2D.Mesh()
MyMesh.generateMeshFromDomain(MyDomain, 0.01)
MyMesh.show()
MyMesh.saveToJson("HughesRuGrandmont")


#Construction de la donnée initiale ---------------------------------------
MyMap = Mesh2D.CellValueMap(MyMesh)
MyMap.generateRandom()
MyMap.show()

opt = dict(model = "hughes",
            filename = "data/HughesRuGrandmont",
            save = True,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001},
            eikoSolver = {  'constrained' : True,
                            'NarrowBandDepth' : 2})

Solver = Splitting.PedestrianSolver(MyMesh, 0.01,0.01, initialDensity = MyMap, options=opt)
#for i in range(5):
#Solver.computeStepsAndShow(2)
#Solver.directions[-1].checkGradientValidity()
#Solver.directions[-1].showVectorField()
Solver.computeUntilEmpty(5000)
