from pathlib import Path

from hughes2d import CellValueMap, EikoSolver, Mesh

MyMesh = Mesh()

#Construction of the mesh
#from hughes2d import NonConvexDomain
#MyDomain = Mesh2D.NonConvexDomain([[0,0],[3,0],[3,3],[0,3]])
#MyDomain.add_exits([[[1,0],[2,0]],[[3,0.5],[3,1.5]],[[0.5,3],[2,3]]])
#MyMesh.generate_mesh_from_domain(MyDomain, 0.01)
#MyMesh.save_to_json("examples/04-Simulating_Eikonal_equation/data/square")

#Loading the mesh
MyMesh.load_from_json(Path(__file__).parent / "data" / "square_mesh.json")

MyMap = CellValueMap(MyMesh)

opt=dict(method = "FMT", constrained = True, NarrowBandDepth = 2)

EikoSolv = EikoSolver(MyMesh, density_map = MyMap, opt = opt)
EikoSolv.compute_field()

#Show as a colored scatter plot
EikoSolv.field_values.show(grid=True, colorscale_name = "magenta")

#Show as a 3D plot
EikoSolv.field_values.show_3d()

#Displays the vector field associated with the computed approximation
EikoSolv.field_values.show_vector_field()
