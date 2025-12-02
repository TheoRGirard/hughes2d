from pathlib import Path

from hughes2d import Mesh, NonConvexDomain

big_square_points = [[0,0],[5,0],[5,5],[0,5]]

small_square_points = [[2,2],[2,3],[3,3],[3,2]]

Exit1 = [[0,1],[0,2]]
Exit2 = [[2.25,3],[2.75,3]]
Exit3 = [[5,1],[5,2]]

Innerwall = [[1,1],[1,4]]

#Domain construction -----------------------------
Domain1 = NonConvexDomain(big_square_points)
Domain1.add_wall(small_square_points, cycle=True)
Domain1.add_wall(Innerwall, cycle=False)
Domain1.add_exit(Exit1)
Domain1.add_exit(Exit2)
Domain1.add_exit(Exit3)

#Definition of the mesh and generation from the domain -------------
MyMesh = Mesh()
MyMesh.generate_mesh_from_domain(Domain1, 0.1)

#Saving the mesh
MyMesh.save_to_json(str(Path(__file__).parent / "square"))

#Loading the mesh
MyMesh.load_from_json(str(Path(__file__).parent / "square_mesh.json"))
