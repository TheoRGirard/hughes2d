from hughes2d import *
import random as alea
import plotly
import matplotlib.pyplot as plt
import numpy as np


Exits = [[[10,0],[10,5]]]

ListOuterWall2 = [[0,0],[10,0],[10,5],[0,5]]


MyDomain = Mesh2D.NonConvexDomain(ListOuterWall2)
MyDomain.addExits(Exits)


dt = 0.02
dx = 0.02
numStep = 1000

lineWidth = 0.2

MyMesh = Mesh2D.Mesh()

#MyMesh.generateMeshFromDomain(MyDomain,dx)
#MyMesh.show()
#MyMesh.saveToJson("data/TestRegu_1")
MyMesh.loadFromJson("data/TestRegu_1_mesh.json")

MyMap = Mesh2D.CellValueMap(MyMesh)
#MyMap.generateRandom()
for i in range(len(MyMap)):
    if(MyMesh.barycenters[i][0] < 9):
        if( int((MyMesh.barycenters[i][1])//lineWidth) % 4 != 0 and int((MyMesh.barycenters[i][0])//lineWidth) % 8 > 1):
            MyMap[i] = 0.95
        else:
            MyMap[i] = 0
    else:
        MyMap[i] = 0

MyMap.show()


opt = dict(model = "hughes",
            filename = "data/TestRegu_1",
            save = False,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.000001},
            eikoSolver = {  'constrained' : True,
                            'NarrowBandDepth' : 2},
            additional_computations = { 'total_mass' : True })

Solver = Splitting.PedestrianSolver(MyMesh, dt,dx, initialDensity = MyMap, costFunction =(lambda x : 1+x), options=opt)

MyVectors = [[1,0] for _ in MyMesh.triangles]

opt2 = dict(model = "constantDirectionField",
            filename = "data/HughesCouloirTmap",
            save = False,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001}
                            )

Solver2 = Splitting.PedestrianSolver(MyMesh, dt,dx, initialDensity = MyMap, directions = MyVectors, options=opt2)

maxDV = [[],[]]
meanDV = [[],[]]
maxDT = [[],[]]
meanDT = [[],[]]

for i in range(numStep):
    DerivV = []
    DerivT = []
    cDerivV = []
    cDerivT = []
    for index,T in enumerate(MyMesh.trianglesWithEdges):
        maxGrad = 0
        maxGradT = 0
        cmaxGrad = 0
        cmaxGradT = 0

        for j,edge1 in enumerate(T):
            for edge2 in T[(j+1):]:
                result = 0
                resultT = 0
                cresult = 0
                cresultT = 0
                for neighbourTriangle in MyMesh.pairsOfTriangles[edge1]:
                    if neighbourTriangle != index:
                        result += ( (Solver.LWRsolver.densityt0[neighbourTriangle] - Solver.LWRsolver.densityt0[index])
                                        /( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])**2
                                            + (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])**2 )
                                        *( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])
                                            *(Solver.directions[index][0])
                                            +(MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])
                                            *(Solver.directions[index][1]))
                                        )
                        resultT += ( (Solver.LWRsolver.densityt0[neighbourTriangle] - Solver.LWRsolver.densityt0[index])
                                        /( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])**2
                                            + (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])**2 )
                                        *( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])
                                            *(Solver.directions[index][1])*(-1)
                                            +(MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])
                                            *(Solver.directions[index][0]))
                                        )
                        cresult += ( (Solver2.LWRsolver.densityt0[neighbourTriangle] - Solver2.LWRsolver.densityt0[index])
                                        /( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])**2
                                            + (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])**2 )
                                        *( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])
                                            *(Solver2.directions[index][0])
                                            +(MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])
                                            *(Solver2.directions[index][1]))
                                        )
                        cresultT += ( (Solver2.LWRsolver.densityt0[neighbourTriangle] - Solver2.LWRsolver.densityt0[index])
                                        /( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])**2
                                            + (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])**2 )
                                        *( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])
                                            *(Solver2.directions[index][1])*(-1)
                                            +(MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])
                                            *(Solver2.directions[index][0]))
                                        )

                for neighbourTriangle in MyMesh.pairsOfTriangles[edge2]:
                    if neighbourTriangle != index:
                        result += ( (Solver.LWRsolver.densityt0[neighbourTriangle] - Solver.LWRsolver.densityt0[index])
                                        /( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])**2
                                            + (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])**2 )
                                        *( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])
                                            *(Solver.directions[index][0])
                                            +(MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])
                                            *(Solver.directions[index][1]))
                                        )
                        resultT += ( (Solver.LWRsolver.densityt0[neighbourTriangle] - Solver.LWRsolver.densityt0[index])
                                        /( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])**2
                                            + (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])**2 )
                                        *( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])
                                            *(Solver.directions[index][1])*(-1)
                                            +(MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])
                                            *(Solver.directions[index][0]))
                                        )
                        cresult += ( (Solver2.LWRsolver.densityt0[neighbourTriangle] - Solver2.LWRsolver.densityt0[index])
                                        /( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])**2
                                            + (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])**2 )
                                        *( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])
                                            *(Solver2.directions[index][0])
                                            +(MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])
                                            *(Solver2.directions[index][1]))
                                        )
                        cresultT += ( (Solver2.LWRsolver.densityt0[neighbourTriangle] - Solver2.LWRsolver.densityt0[index])
                                        /( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])**2
                                            + (MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])**2 )
                                        *( (MyMesh.barycenters[neighbourTriangle][0] - MyMesh.barycenters[index][0])
                                            *(Solver2.directions[index][1])*(-1)
                                            +(MyMesh.barycenters[neighbourTriangle][1] - MyMesh.barycenters[index][1])
                                            *(Solver2.directions[index][0]))
                                        )
                if abs(result) > maxGrad:
                    maxGrad = abs(result)
                if abs(resultT) > maxGradT:
                    maxGradT = abs(resultT)
                if abs(cresult) > cmaxGrad:
                    cmaxGrad = abs(cresult)
                if abs(cresultT) > cmaxGradT:
                    cmaxGradT = abs(cresultT)

        DerivV.append(maxGrad)
        DerivT.append(maxGradT)
        cDerivV.append(cmaxGrad)
        cDerivT.append(cmaxGradT)

    maxDV[0].append(max(DerivV))
    meanDV[0].append(sum(DerivV)/len(DerivV))
    maxDV[1].append(max(cDerivV))
    meanDV[1].append(sum(cDerivV)/len(cDerivV))
    maxDT[0].append(max(DerivT))
    meanDT[0].append(sum(DerivT)/len(DerivT))
    maxDT[1].append(max(cDerivT))
    meanDT[1].append(sum(cDerivT)/len(cDerivT))
    Solver.computeStep()
    Solver2.computeStep()
    print("Num step : ", i ,"/", numStep-1)

PlotDatas = [[maxDV,maxDT],[meanDV,meanDT]]

T = [i*dt for i in range(numStep)]
models =["direction of V", "orthogonal to V"]
for j,plotName in enumerate(["Max norm", "Mean norm"]):
    fig, axs = plt.subplots(1)
    for i,model in enumerate(models):
        axs.plot(T,PlotDatas[j][i][0],"C"+str(i)+"-", label = "Hughes with "+ model)
        axs.plot(T,PlotDatas[j][i][1],"C"+str(i)+"--", label = "LWR constant field with "+model)
    axs.set_xlabel("t")
    axs.set_yscale("log")
    axs.set_xlim(0,25)
    axs.set_ylim(0.1,10)
    axs.legend()
    #plt.show()
    plt.savefig("figs/Compare5Log"+plotName[:3]+str(dx)+".png")

for j,plotName in enumerate(["Max norm", "Mean norm"]):
    fig, axs = plt.subplots(1)
    for i,model in enumerate(models):
        axs.plot(T,PlotDatas[j][i][0],"C"+str(i)+"-", label = "Hughes with "+ model)
        axs.plot(T,PlotDatas[j][i][1],"C"+str(i)+"--", label = "LWR constant field with "+model)
    axs.set_xlabel("t")
    axs.set_xlim(0,25)
    axs.set_ylim(0,10)
    axs.legend()
    #plt.show()
    plt.savefig("figs/Compare5"+plotName[:3]+str(dx)+".png")
