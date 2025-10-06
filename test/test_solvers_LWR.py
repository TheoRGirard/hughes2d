from hughes2d.Mesh2D import *
from hughes2d.LWR2D import *
import pytest

"""
This file tests the LWR solver.
"""

MyDomain = NonConvexDomain([[0,0],[10,0],[10,5],[0,5]])
MyDomain.addExit([[10,0],[10,5]])

MyMesh = Mesh()

#Mesh generation
MyMesh.generateMeshFromDomain(MyDomain, 0.2)
print("Number of triangles : ", len(MyMesh.triangles))
MyMesh.saveToJson("test/ressources/test_LWR")

#Mesh loading
#MyMesh.loadFromJson("test/ressources/test_LWR_mesh.json")

VectorField = [[1,0] for _ in MyMesh.triangles]

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

@pytest.mark.parametrize("method, convexFlux, dt, CFL, numSteps, precision", [("midVector",True, 0.3, False, 0, 0.005),("midVector",True, 0.1, False, 100, 0.05),("midVector",False, 0.02, True, 100, 0.03),("tmap", True, 0.1, False, 100, 0.05),("tmap", False, 0.02, True, 100, 0.03)])
def test_LWR(method, convexFlux, dt, CFL, numSteps, precision):
    options = dict(method = method, convexFlux = convexFlux, debugging = True)
    if not convexFlux:
        options['ApproximationThreshold'] = 0.001
    MySolver = LWRSolver(MyMesh, dt = dt, previousDensity = InitialDatum, directionMap = VectorField, opt=options)

    assert MySolver.checkCFL() == CFL

    for i in range(numSteps):
        slice = CellValueMap(MyMesh)
        slice.values = [explicitSolFunc((i+1)*dt, P[0]) for P in MyMesh.barycenters]

        MySolver.computeNextStep()

        meanDiff = 0
        for j in range(len(MyMesh.triangles)):
            meanDiff += abs(MySolver.densityt1[j] - slice[j])*MyMesh.cellAreas[j]
        meanDiff = meanDiff/50

        assert meanDiff < precision

        MySolver.update(VectorField)
