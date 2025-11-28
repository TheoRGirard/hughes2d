
"""
This example shows how to compare different simulation. Here we compare the
approximation of Hughes model with the explicit solution.
We also compare the constant direction field simulation when the constant vector field
is set as the explicit vector field.
It is interesting to notice that the Hughes' scheme is closer to the explicit solution
than the LWR + explicit vectors simulation.
"""
from pathlib import Path

import random as alea

import matplotlib.pyplot as plt
import numpy as np
import plotly

from hughes2d import CellValueMap, Mesh, NonConvexDomain, PedestrianSolver

#Parameters of the simulation
dt = 0.025
dx = 0.05
num_step = 2000

#Construction of the mesh
Exits = [[[10,0],[10,5]]]
ListOuterWall2 = [[0,0],[10,0],[10,5],[0,5]]
MyDomain = NonConvexDomain(ListOuterWall2)
MyDomain.add_exits(Exits)
MyMesh = Mesh()
MyMesh.generate_mesh_from_domain(MyDomain,dx)

#Construction of the initial datum
initial_datum = CellValueMap(MyMesh)
for i in range(len(initial_datum)):
    if(MyMesh.barycenters[i][0] < 5):
        initial_datum[i] = 0.7
    else:
        initial_datum[i] = 0


#Computation of the explicit solution ----------------------------
def explicit_sol_func(t,x):
    if (x <= 0.3*t and t <= 5/0.7) or (x <= 5 + t - 2 *np.sqrt(5*0.7*t) and t > 5/0.7):
        return 0
    return min(0.7,max(0,(5-x)/(2*t) + 0.5))

ExplicitSol = [initial_datum.values]
for i in range(1,num_step):
    actual_slice = CellValueMap(MyMesh)
    actual_slice.values = [explicit_sol_func(i*dt, P[0]) for P in MyMesh.barycenters]
    ExplicitSol.append(actual_slice)

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

SolverConstantDirection = PedestrianSolver(MyMesh, dt,
                                                     initial_density = initial_datum,
                                                     directions = ExplicitVectorField,
                                                     options=opt)


#Constructing the hughes solver
opt2 = dict(model = "hughes",
            save = False,
            verbose = True,
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001}
                            )

SolverHughes = PedestrianSolver(MyMesh, dt,
                                          initial_density = initial_datum,
                                          options=opt2)


#Computing the L1 differences with the explicit solution step by step
L1Diffs = [[],[]]
for j in range(num_step):
    SolverConstantDirection.computeStep()
    SolverHughes.computeStep()

    D = 0
    D2 = 0
    for k in range(len(MyMesh.triangles)):
        D += (abs(SolverConstantDirection.lwr_solver.densityt1[k]-ExplicitSol[j][k])
              *MyMesh.cell_areas[k])
        D2 += (abs(SolverHughes.lwr_solver.densityt1[k]-ExplicitSol[j][k])
               *MyMesh.cell_areas[k])
    L1Diffs[0].append(D/50)
    L1Diffs[1].append(D2/50)

#Plotting both graphs
T = [i*dt for i in range(num_step)]
models =["constant explicit direction field", "hughes"]
fig, axs = plt.subplots(1)
for i,model in enumerate(models):
    axs.plot(T,L1Diffs[i], label = model)
axs.set_title("L¹ difference with the explicit solution with |Δ| = "+str(dx)+".")
axs.set_xlabel("t")
axs.legend()

plt.savefig(Path(__file__).parent / "figs" / "CompareExplicit.png")
