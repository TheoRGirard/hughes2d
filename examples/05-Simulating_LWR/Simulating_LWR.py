from hughes2d import CellValueMap, LWRSolver, Mesh, NonConvexDomain

#Parameters of the simulation
dt = 0.025
dx = 0.05
num_step = 200

#Construction of the mesh
MyDomain = NonConvexDomain([[0,0],[10,0],[10,5],[0,5]])
MyDomain.add_exits([[[10,0],[10,5]]])
MyMesh = Mesh()
MyMesh.generate_mesh_from_domain(MyDomain,dx)

#Construction of an initial datum
InitialDatum = CellValueMap(MyMesh)
InitialDatum.generate_random()

#Plot the InitialDatum
InitialDatum.show(preference="matplotlib")

#Construction of a vector field
VectorField = [[1,0] for _ in MyMesh.triangles]

#Construction of the LWR solver
opt = dict( method = "midVector",
            convexFlux = True )
MySolver = LWRSolver(MyMesh, dt, previous_density=InitialDatum,
                     direction_map =VectorField, opt = opt)

#Compute the number of steps num_step
for _ in range(num_step):
    MySolver.compute_next_step()
    MySolver.update(VectorField)

#Plot the final state of the density
MySolver.show_density(t=1, preference="matplotlib")
