from hughes2d.Mesh2D import *
from hughes2d.Splitting import *
import pytest

try:
    import plotly.figure_factory as ff
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
#MyMesh.saveToJson("test/ressources/test_pedestrian")

#Mesh loading
MyMesh.loadFromJson("test/ressources/test_pedestrian_mesh.json")

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

    MySolver = PedestrianSolver(MyMesh, dt = dt, initialDensity = InitialDatum, options=options)

    for i in range(numSteps):
        slice = CellValueMap(MyMesh)
        slice.values = [explicitSolFunc((i+1)*dt, P[0]) for P in MyMesh.barycenters]

        MySolver.computeStep()

        meanDiff = 0
        for j in range(len(MyMesh.triangles)):
            meanDiff += abs(MySolver.LWRsolver.densityt0[j] - slice[j])*MyMesh.cellAreas[j]
        meanDiff = meanDiff/50

        assert meanDiff < precision

def test_pedestrian_save():
    options = dict(model = "constantDirectionField",
                    save = True,
                    filename = "test/test_pedestrian_outputs/test",
                    additional_computations = dict(total_mass = True, max_density = True)
                    )

    InitialDatum2  = CellValueMap(MyMesh)
    InitialDatum2.setConstantCircle([9.5,2.5],0.3,0.3)

    MySolver = PedestrianSolver(MyMesh, dt = 0.1, initialDensity = InitialDatum2, options=options)

    MySolver.computeUntilEmpty()

    assert MySolver.timeStep < 20

    MySolver.saveToJson()

def test_pedestrian_plot():
    if(go):
        options = dict(model = "hughes")

        MySolver = PedestrianSolver(MyMesh, dt = 0.1, initialDensity = InitialDatum, options=options)

        MySolver.computeStepsAndShow(10)
