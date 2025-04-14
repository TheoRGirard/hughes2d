from hughes2d import *
import random as alea
import plotly
import matplotlib.pyplot as plt
import numpy as np


Exits = [[[10,0],[10,5]]]

ListOuterWall2 = [[0,0],[10,0],[10,5],[0,5]]


MyDomain = Mesh2D.NonConvexDomain(ListOuterWall2)
MyDomain.addExits(Exits)


dt = 0.025
dx = 0.05
numStep = 2000

MyMesh = Mesh2D.Mesh()

MyMesh.generateMeshFromDomain(MyDomain,dx)
#MyMesh.show()


MyMap = Mesh2D.CellValueMap(MyMesh)
#MyMap.generateRandom()
for i in range(len(MyMap)):
    if(MyMesh.barycenters[i][0] < 5):
        MyMap[i] = 0.7
    else:
        MyMap[i] = 0



def explicitSolFunc(t,x):
    if x <= 0.3*t and t <= 5/0.7:
        return 0
    elif x <= 5 + t - 2 *np.sqrt(5*0.7*t) and t > 5/0.7:
        return 0
    return min(0.7,max(0,(5-x)/(2*t) + 0.5))

ExplicitSol = []
for i in range(numStep):
    slice = Mesh2D.CellValueMap(MyMesh)
    slice.values = [explicitSolFunc(i*dt, P[0]) for P in MyMesh.barycenters]
    ExplicitSol.append(slice)

MyEmptyMap = Mesh2D.CellValueMap(MyMesh)

Opt3=dict(method = "FMT", constrained = True, NarrowBandDepth = 2)
FMTCSolv = EikonalSolver.EikoSolver(MyMesh,DensityMap = MyEmptyMap,opt=Opt3)
FMTCSolv.computeField()
MyVectors2 = FMTCSolv.fieldValues.computeGradientFlow()

MyVectors = [[1,0] for _ in MyMesh.triangles]

opt = dict(model = "constantDirectionField",
            filename = "data/HughesCouloirTmap",
            save = False,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001}
                            )

Solver = Splitting.PedestrianSolver(MyMesh, dt,dx, initialDensity = MyMap, directions = MyVectors, options=opt)

opt3 = dict(model = "constantDirectionField",
            filename = "data/HughesCouloirTmap",
            save = False,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001}
                            )

Solver3 = Splitting.PedestrianSolver(MyMesh, dt,dx, initialDensity = MyMap, directions = MyVectors2, options=opt3)

opt2 = dict(model = "hughes",
            filename = "data/HughesCouloirTmap",
            save = False,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001},
            eikoSolver = { 'method' : 'FMT',
                            'constrained' : True,
                            'NarrowBandDepth' : 2
                            })

Solver2 = Splitting.PedestrianSolver(MyMesh, dt,dx, initialDensity = MyMap, options=opt2)

L1Diffs = [[],[],[]]
for j in range(numStep):
    Solver.computeStep()
    Solver2.computeStep()
    Solver3.computeStep()
    if(j == 800):
        Map = Mesh2D.CellValueMap(MyMesh)
        Map.values = ExplicitSol[j]
        Map.show()
        Map.values = Solver.LWRsolver.densityt1
        Map.show()
        Map.values = [np.abs(Solver.LWRsolver.densityt1[k]-ExplicitSol[j][k]) for k in range(len(MyMesh.triangles))]
        Map.show()
        print(max(Map.values))
        Map.values = Solver2.LWRsolver.densityt1
        Map.show()
        Map.values = [np.abs(Solver2.LWRsolver.densityt1[k]-ExplicitSol[j][k]) for k in range(len(MyMesh.triangles))]
        Map.show()
        print(max(Map.values))


    D = 0
    D2 = 0
    D3 = 0
    for k in range(len(MyMesh.triangles)):
        D += abs(Solver.LWRsolver.densityt1[k]-ExplicitSol[j][k])*MyMesh.cellAreas[k]
        D2 += abs(Solver2.LWRsolver.densityt1[k]-ExplicitSol[j][k])*MyMesh.cellAreas[k]
        D3 += abs(Solver3.LWRsolver.densityt1[k]-ExplicitSol[j][k])*MyMesh.cellAreas[k]
    L1Diffs[0].append(D/50)
    L1Diffs[1].append(D2/50)
    L1Diffs[2].append(D3/50)

T = [i*dt for i in range(numStep)]
models =["Hughes' model", "LWR - explicit vectors", "LWR - FMTC vectors"]
fig, axs = plt.subplots(1)
for i,model in enumerate(models):
    axs.plot(T,L1Diffs[i], label = model)
axs.set_title("L¹ difference with the explicit solution with |Δ| = "+str(dx)+".")
axs.set_xlabel("t")
axs.legend()
#plt.show()
plt.savefig("figs/CompareExplicit.png")
