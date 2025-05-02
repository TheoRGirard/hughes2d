from hughes2d import *
import random as alea
import plotly
import matplotlib.pyplot as plt
import numpy as np


Exits = [[[10,0],[10,5]]]

ListOuterWall2 = [[0,0],[10,0],[10,5],[0,5]]


MyDomain = Mesh2D.NonConvexDomain(ListOuterWall2)
MyDomain.addExits(Exits)


dt = 0.04
dx = 0.04
numStep = 1000

lineWidth = 0.5

MyMesh = Mesh2D.Mesh()

MyMesh.generateMeshFromDomain(MyDomain,dx)
#MyMesh.show()


MyMap = Mesh2D.CellValueMap(MyMesh)
#MyMap.generateRandom()
for i in range(len(MyMap)):
    if(MyMesh.barycenters[i][0] < 9):
        if( int(MyMesh.barycenters[i][1]//lineWidth) % 2 == 0 and int(MyMesh.barycenters[i][0]//lineWidth) % 4 != 0):
            MyMap[i] = 0.9
        else:
            MyMap[i] = 0
    else:
        MyMap[i] = 0

MyMap.show()


opt = dict(model = "hughes",
            filename = "data/TestRegu_"+str(np.round(dx,1)),
            save = False,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.000001},
            eikoSolver = {  'constrained' : True,
                            'NarrowBandDepth' : 2},
            additional_computations = { 'total_mass' : True })

Solver = Splitting.PedestrianSolver(MyMesh, 0.015,0.03, initialDensity = MyMap, costFunction =(lambda x : 1+x), options=opt)

maxDV = []
meanDV = []
maxDT = []
meanDT = []

for i in range(numStep):
    DerivV = []
    DerivT = []
    for index,T in enumerate(MyMesh.trianglesWithEdges):
        maxGrad = 0
        maxGradT = 0
        Visited = [index]
        for j,edge1 in enumerate(T):
            for edge2 in T[(j+1):]:
                result = 0
                resultT = 0
                for neighbourTriangle in MyMesh.pairsOfTriangles[edge1]:
                    if neighbourTriangle != index:
                        result += ( (Solver.LWRsolver.densityt0[neighbourTriangle] - Solver.LWRsolver.densityt0[index])
                                        /( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])**2
                                            + (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])**2 )
                                        *( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0]) )
                                        )
                        resultT += ( (Solver.LWRsolver.densityt0[neighbourTriangle] - Solver.LWRsolver.densityt0[index])
                                        /( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])**2
                                            + (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])**2 )
                                        *( (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1]) )
                                        )

                for neighbourTriangle in MyMesh.pairsOfTriangles[edge2]:
                    if neighbourTriangle != index:
                        result += ( (Solver.LWRsolver.densityt0[neighbourTriangle] - Solver.LWRsolver.densityt0[index])
                                        /( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])**2
                                            + (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])**2 )
                                        *( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0]) )
                                        )
                        resultT += ( (Solver.LWRsolver.densityt0[neighbourTriangle] - Solver.LWRsolver.densityt0[index])
                                        /( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])**2
                                            + (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])**2 )
                                        *( (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1]) )
                                        )
                if abs(result) > maxGrad:
                    maxGrad = abs(result)
                if abs(resultT) > maxGradT:
                    maxGradT = abs(resultT)

        DerivV.append(maxGrad)
        DerivT.append(maxGradT)

    maxDV.append(max(DerivV))
    meanDV.append(sum(DerivV)/len(DerivV))
    maxDT.append(max(DerivT))
    meanDT.append(sum(DerivT)/len(DerivT))
    Solver.computeStep()
    print("Num step : ", i ,"/", numStep-1)

PlotDatas = [[maxDV,maxDT],[meanDV,meanDT]]

T = [i*dt for i in range(numStep)]
models =["direction V", "orthogonal to V"]
for j,plotName in enumerate(["Max norm", "Mean norm"]):
    fig, axs = plt.subplots(1)
    for i,model in enumerate(models):
        axs.plot(T,PlotDatas[j][i], label = model)
    axs.set_title(plotName + " of the spatial derivative of the density depending on the direction.")
    axs.set_xlabel("t")
    axs.set_yscale("log")
    axs.set_ylimits(0.01,10)
    axs.legend()
    #plt.show()
    plt.savefig("figs/CompareDistance2"+plotName[:3]+str(dx)+".png")
