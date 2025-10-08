
"""
This example shows how to compare different simulation. Here we compare the approximation of Hughes model with the explicit solution.
We also compare the constant direction field simulation when the constant vector field is set as the explicit vector field.
It is interesting to notice that the Hughes' scheme is closer to the explicit solution than the LWR + explicit vectors simulation.
"""


from hughes2d import *
import random as alea
import plotly
import matplotlib.pyplot as plt
import numpy as np

#Parameters of the simulation
dt = 0.025
dx = 0.05
numStep = 2000

#Construction of the mesh
Exits = [[[10,0],[10,5]]]
ListOuterWall2 = [[0,0],[10,0],[10,5],[0,5]]
MyDomain = Mesh2D.NonConvexDomain(ListOuterWall2)
MyDomain.addExits(Exits)
MyMesh = Mesh2D.Mesh()
MyMesh.generateMeshFromDomain(MyDomain,dx)

#Construction of the initial datum
InitialDatum = Mesh2D.CellValueMap(MyMesh)
for i in range(len(InitialDatum)):
    if(MyMesh.barycenters[i][0] < 5):
        InitialDatum[i] = 0.7
    else:
        InitialDatum[i] = 0


#Computation of the explicit solution ----------------------------
def explicitSolFunc(t,x):
    if x <= 0.3*t and t <= 5/0.7:
        return 0
    elif x <= 5 + t - 2 *np.sqrt(5*0.7*t) and t > 5/0.7:
        return 0
    return min(0.7,max(0,(5-x)/(2*t) + 0.5))

ExplicitSol = [InitialDatum.values]
for i in range(1,numStep):
    slice = Mesh2D.CellValueMap(MyMesh)
    slice.values = [explicitSolFunc(i*dt, P[0]) for P in MyMesh.barycenters]
    ExplicitSol.append(slice)

#Construction of the explicit vectors
ExplicitVectorField = [[1,0] for _ in MyMesh.triangles]

#Constructing the constant direction field solver
opt = dict(model = "constantDirectionField",
            save = False,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001}
                            )

SolverConstantDirection = Splitting.PedestrianSolver(MyMesh, dt, initialDensity = InitialDatum, directions = ExplicitVectorField, options=opt)


#Constructing the hughes solver
opt2 = dict(model = "hughes",
            save = False,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001}
                            )

SolverHughes = Splitting.PedestrianSolver(MyMesh, dt, initialDensity = InitialDatum, options=opt2)


#Computing the L1 differences with the explicit solution step by step
L1Diffs = [[],[]]
for j in range(numStep):
    SolverConstantDirection.computeStep()
    SolverHughes.computeStep()

    D = 0
    D2 = 0
    for k in range(len(MyMesh.triangles)):
        D += abs(SolverConstantDirection.LWRsolver.densityt1[k]-ExplicitSol[j][k])*MyMesh.cellAreas[k]
        D2 += abs(SolverHughes.LWRsolver.densityt1[k]-ExplicitSol[j][k])*MyMesh.cellAreas[k]
    L1Diffs[0].append(D/50)
    L1Diffs[1].append(D2/50)

#Plotting both graphs
T = [i*dt for i in range(numStep)]
models =["constant explicit direction field", "hughes"]
fig, axs = plt.subplots(1)
for i,model in enumerate(models):
    axs.plot(T,L1Diffs[i], label = model)
axs.set_title("L¹ difference with the explicit solution with |Δ| = "+str(dx)+".")
axs.set_xlabel("t")
axs.legend()

plt.savefig("examples/03-Simulating_Hughes_2/figs/CompareExplicit.png")
