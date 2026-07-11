from hughes2d import Plotter
from hughes2d.EikonalSolver import EikoSolver
from hughes2d.LWR2D import LWRSolver
from hughes2d.Mesh2D import CellValueMap, Mesh, NonConvexDomain, VertexValueMap
from hughes2d.Splitting import PedestrianSolver

__all__ = [
        "CellValueMap",
        "EikoSolver",
        "LWRSolver",
        "Mesh",
        "NonConvexDomain",
        "PedestrianSolver",
        "Plotter",
        "VertexValueMap",
        ]

__version__ = '1.1.6'
