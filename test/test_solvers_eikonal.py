import pytest

from hughes2d.EikonalSolver import EikoSolver
from hughes2d.Mesh2D import CellValueMap, Mesh

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
MyMesh.load_from_json("test/ressources/test_Eiko_mesh.json")

InitialDatum = CellValueMap(MyMesh)
InitialDatum.set_constant(0)



@pytest.mark.parametrize("method,NBdepth,constrained, precision", [("FMT",2,True, 0.018),("FMT",1,True, 0.4),("FMT",2,False, 0.07),("FMT",1,False, 0.005),("FME",1,False, 0.06)])
def test_simple_solution_ID0(method,NBdepth,constrained, precision):
    InitialDatum.set_constant(0)
    options = dict(method = method, NarrowBandDepth = NBdepth, constrained = constrained)
    MySolver = EikoSolver(MyMesh, density_map = InitialDatum, opt=options)
    MySolver.compute_field()

    maxError = 0
    for i,point in enumerate(MyMesh.vertices):
        error = abs( (1 - point[0]) - MySolver.field_values[i] )
        if error > maxError :
            maxError = error

    assert maxError < precision

@pytest.mark.parametrize("method,NBdepth,constrained, precision", [("FMT",2,True, 0.06),("FMT",1,True, 1),("FMT",2,False, 0.25),("FMT",1,False, 0.005),("FME",1,False, 0.18)])
def test_simple_solution_ID2(method,NBdepth,constrained, precision):
    InitialDatum.set_constant(2)
    options = dict(method = method, NarrowBandDepth = NBdepth, constrained = constrained)
    MySolver = EikoSolver(MyMesh, density_map = InitialDatum, opt=options)
    MySolver.compute_field()

    maxError = 0
    for i,point in enumerate(MyMesh.vertices):
        error = abs( (3 - 3*point[0]) - MySolver.field_values[i] )
        if error > maxError :
            maxError = error

    assert maxError < precision

def test_plot_narrowBand():
    InitialDatum.set_constant(1)
    options = dict(method = "FME")
    MySolver = EikoSolver(MyMesh, density_map = InitialDatum, opt=options)
    if go:
        MySolver.show_narrow_band_after_step(10)
    else:
        with pytest.raises(ImportError):
            MySolver.show_narrow_band_after_step(10)
