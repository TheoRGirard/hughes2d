from hughes2d.Mesh2D import *
from hughes2d.EikonalSolver import *
import pytest

"""
This file tests the eikonal solver.
"""

try:
    import plotly.figure_factory as ff
    import plotly.graph_objects as go
except ImportError:
    go = None

MyMesh = Mesh()

#Mesh generation
#MyDomain = NonConvexDomain([[0,0],[0,1],[1,1],[1,0]])
#MyDomain.addExit([[1,0],[1,1]])
#MyMesh.generateMeshFromDomain(MyDomain, 0.005)
#MyMesh.saveToJson("test/ressources/test_Eiko")

#Mesh loading
MyMesh.loadFromJson("test/ressources/test_Eiko_mesh.json")

InitialDatum = CellValueMap(MyMesh)
InitialDatum.setConstant(0)



@pytest.mark.parametrize("method,NBdepth,constrained, precision", [("FMT",2,True, 0.018),("FMT",1,True, 0.4),("FMT",2,False, 0.07),("FMT",1,False, 0.005),("FME",1,False, 0.06)])
def test_simple_solution_ID0(method,NBdepth,constrained, precision):
    InitialDatum.setConstant(0)
    options = dict(method = method, NarrowBandDepth = NBdepth, constrained = constrained)
    MySolver = EikoSolver(MyMesh, DensityMap = InitialDatum, opt=options)
    MySolver.computeField()

    maxError = 0
    for i,point in enumerate(MyMesh.vertices):
        error = abs( (1 - point[0]) - MySolver.fieldValues[i] )
        if error > maxError :
            maxError = error

    assert maxError < precision

@pytest.mark.parametrize("method,NBdepth,constrained, precision", [("FMT",2,True, 0.06),("FMT",1,True, 1),("FMT",2,False, 0.25),("FMT",1,False, 0.005),("FME",1,False, 0.18)])
def test_simple_solution_ID2(method,NBdepth,constrained, precision):
    InitialDatum.setConstant(2)
    options = dict(method = method, NarrowBandDepth = NBdepth, constrained = constrained)
    MySolver = EikoSolver(MyMesh, DensityMap = InitialDatum, opt=options)
    MySolver.computeField()

    maxError = 0
    for i,point in enumerate(MyMesh.vertices):
        error = abs( (3 - 3*point[0]) - MySolver.fieldValues[i] )
        if error > maxError :
            maxError = error

    assert maxError < precision

def test_plot_narrowBand():
    InitialDatum.setConstant(1)
    options = dict(method = "FME")
    MySolver = EikoSolver(MyMesh, DensityMap = InitialDatum, opt=options)
    if go:
        MySolver.showNarrowBandAfterStep(10)
    else:
        with pytest.raises(ImportError):
            MySolver.showNarrowBandAfterStep(10)
