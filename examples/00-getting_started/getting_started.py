from pathlib import Path

from hughes2d import CellValueMap, Mesh, NonConvexDomain, PedestrianSolver

file_path = str(Path(__file__).parent / "gettingStartedSimu")

#Construction of the domain--------------------------------
MyDomain = NonConvexDomain([[0,0],[0,1],[1,1],[1,0]])
MyDomain.add_exits([[[1,0],[1,1]]])

#construction of the Mesh--------------------------------------
MyMesh = Mesh()
MyMesh.generate_mesh_from_domain(MyDomain, 0.01)
MyMesh.save_to_json(file_path)


#Construction of the initial datum ---------------------------------------
MyMap = CellValueMap(MyMesh)
MyMap.generate_random()

#Setting the options for the simulation-----------------------------------------
opt = {
    "model" : "hughes",
    "filename" : file_path,
    "save" : True,
    "verbose" : False,
        }

#Creating the solver and computing---------------------------------------------------
Solver = PedestrianSolver(MyMesh, 0.01, initial_density = MyMap, options=opt)
Solver.compute_until_empty(100)

#Converting the data to a mp4 video------------------------------------------
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

if plt:
    from hughes2d import Plotter
Plotter.convert_to_mp4(file_path, limits=[[0,1],[0,1]])
