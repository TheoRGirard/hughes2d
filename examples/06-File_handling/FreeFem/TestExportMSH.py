from hughes2d import Mesh, NonConvexDomain

bigSquarePoints = [[0,0],[5,0],[5,5],[0,5]]
smallSquarePoints = [[2,2],[2,3],[3,3],[3,2]]
Exit1 = [[0,1],[0,2]]
Exit2 = [[2.25,3],[2.75,3]]
InWall = [[1,1],[3,3]]

Domain1 = NonConvexDomain(bigSquarePoints)
Domain1.add_wall(smallSquarePoints, cycle=True)
Domain1.add_exit(Exit1)
Domain1.add_exit(Exit2)

Mesh1 = Mesh()
Mesh1.generate_mesh_from_domain(Domain1, 0.1)
Mesh1.show()

Mesh1.export_mesh_msh_free_fem("examples/06-File_handling/FreeFem/simple_domain_FF.msh")
