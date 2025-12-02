import numpy as np
import pytest

from hughes2d.Mesh2D import CellValueMap, Mesh
from hughes2d.Splitting import PedestrianSolver

try:
    import plotly.graph_objects as go
except ImportError:
    go = None

"""
This file tests the LWR solver.
"""

MyMesh = Mesh()

#Mesh generation
#MyDomain = NonConvexDomain([[0,0],[10,0],[10,5],[0,5]])
#MyDomain.addExit([[10,0],[10,5]])
#MyMesh.generateMeshFromDomain(MyDomain, 0.2)
#MyMesh.save_to_json("test/ressources/test_pedestrian")

#Mesh loading
MyMesh.load_from_json("test/ressources/test_pedestrian_mesh.json")

InitialDatum  = CellValueMap(MyMesh)
for i in range(len(InitialDatum)):
    if(MyMesh.barycenters[i][0] < 5):
        InitialDatum [i] = 0.7
    else:
        InitialDatum [i] = 0


def explicitSolFunc(t,x):
    if x <= 0.3*t and t <= 5/0.7:
        return 0
    elif x <= 5 + t - 2 *np.sqrt(5*0.7*t) and t > 5/0.7:
        return 0
    return min(0.7,max(0,(5-x)/(2*t) + 0.5))

@pytest.mark.parametrize("model, dt, numSteps, precision", [("hughes", 0.1, 100, 0.04),("colombo-garavello", 0.1, 100, 0.09),("constantDirectionField", 0.1, 100, 0.09)])
def test_pedestrian_models(model, dt, numSteps, precision):
    options = dict(model = model)

    MySolver = PedestrianSolver(MyMesh, dt = dt, initial_density = InitialDatum, options=options)

    for i in range(numSteps):
        slice = CellValueMap(MyMesh)
        slice.values = [explicitSolFunc((i+1)*dt, P[0]) for P in MyMesh.barycenters]

        MySolver.compute_step()

        meanDiff = 0
        for j in range(len(MyMesh.triangles)):
            meanDiff += abs(MySolver.lwr_solver.densityt0[j] - slice[j])*MyMesh.cell_areas[j]
        meanDiff = meanDiff/50

        assert meanDiff < precision

def test_pedestrian_save():
    options = dict(model = "constantDirectionField",
                    save = True,
                    filename = "test/test_pedestrian_outputs/test",
                    additional_computations = dict(total_mass = True, max_density = True)
                    )

    InitialDatum2  = CellValueMap(MyMesh)
    InitialDatum2.set_constant_circle([9.5,2.5],0.3,0.3)

    MySolver = PedestrianSolver(MyMesh, dt = 0.1, initial_density = InitialDatum2, options=options)

    MySolver.compute_until_empty()

    assert MySolver.time_step < 20

    MySolver.save_to_json()

def test_pedestrian_plot():
    if(go):
        options = dict(model = "hughes")

        MySolver = PedestrianSolver(MyMesh, dt = 0.1, initial_density = InitialDatum, options=options)

        MySolver.compute_steps_and_show(10)
