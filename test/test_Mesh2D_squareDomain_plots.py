from hughes2d.Mesh2D import *
import pytest

"""
This file tests the plotting methods for the specific example tested in the file "test_Mesh2D_squareDomain.py"
"""

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import matplotlib.cm as cm
    import matplotlib.collections as collections
except ImportError:
    plt = None

try:
    import plotly.figure_factory as ff
    import plotly.graph_objects as go
except ImportError:
    go = None


bigSquarePoints = [[0,0],[5,0],[5,5],[0,5]]
smallSquarePoints = [[2,2],[2,3],[3,3],[3,2]]
Exit1 = [[0,1],[0,2]]
Exit2 = [[2.25,3],[2.75,3]]
InWall = [[1,1],[3,3]]

Domain1 = NonConvexDomain(bigSquarePoints)
Domain1.addWall(smallSquarePoints, cycle=True)
Domain1.addBoundaryPoint([2,5])
Domain1.addWallPoint([2.5,3])
Domain1.addExits([Exit1, Exit2])

Mesh1 = Mesh()
Mesh1.generateMeshFromDomain(Domain1,0.1)

Map1 = CellValueMap(Mesh1)
for i in range(len(Map1)):
    Map1[i] = 0.1

#Plotting method to check the resulting domain:
def test_show_Domain():
    if go:
        Domain1.show(preference="plotly")
    if plt:
        Domain1.show(preference="matplotlib")
    if not go and not plt:
        with pytest.raises(ImportError):
            Domain1.show()

#Plotting method to check the resulting mesh:
def test_show_Mesh():
    if go:
        Mesh1.show(preference="plotly")
    if plt:
        Mesh1.show(preference="matplotlib")
    if not go and not plt:
        with pytest.raises(ImportError):
            Mesh1.show()

def test_show_cellValueMap():
    Map1.setConstantCircle(center=[1,1], radius = 1, value = 0.5)
    if go:
        Map1.show(preference="plotly")
    if plt:
        Map1.show(preference="matplotlib")
    if not go and not plt:
        with pytest.raises(ImportError):
            Map1.show()
