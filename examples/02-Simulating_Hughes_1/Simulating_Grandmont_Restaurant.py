from pathlib import Path

from hughes2d import CellValueMap, Mesh, NonConvexDomain, PedestrianSolver

"""
This example illustrates how to simulate pedestrian flows on a complex domain that
represents the university restaurant of Tours.
"""

Exits = [[[2,3.5],[1.5,3.5]],[[5.25,4.25],[5.75,4.25]],[[9.75,4.25],[9.75,4.75]]]

ListOuterWall2 = [[8.5,10.75],[9,10],[8.25,9.5],[9.25,8.25],[10,8.75],[10.25,8.5],
                    [9.75,7.75],[10.25,7.25],[10.75,7.75],[12.25,5.5],[13,4],[12,3.75],[10.25,3.50],
                    [8.25,3.25],[6.25,3.25],[6.25,4.25],[5.75,4.25],[5.25,4.25],[5.25,3.25],
                    [4.75,3.25],[2,3.5],[1.5,3.5],[1,3.5],[2,5.25],[4,8.25],[5.25,7.25],
                    [6.25,8.5],[6.25,9.5],[5.5,10],[6,10.75],[6.25,10.5],[6.5,10.25],
                    [7,10],[7.5,10],[8,10.25],[8.25,10.5]]

Entries = [[[6.25,10.5],[6.5,10.25]],[[8,10.25],[8.25,10.5]]]

InnerWall1 = [[7,9.5],[7.75,9.5],[8,9.25],[8,8.75],[7.75,8.5],[7,8.5]]

InnerWall2 = [[6.25,8],[8.75,8],[9.25,7.25],[9.5,7.5],[9.75,7.25],[9.25,6.75],
              [10.25,5.5],[10.25,5],[10,4.75],[9.75,4.75],[9.75,4.25],[8.25,4.25],
                [8.25,4.5],[5,4.5],[4.25,4.75],[4.25,5.25],[5.75,7.25]]

#Construction of the domain--------------------------------
MyDomain = NonConvexDomain(ListOuterWall2)
MyDomain.add_wall(InnerWall1,cycle=True)
MyDomain.add_wall(InnerWall2,cycle=True)
MyDomain.add_exits(Exits)
MyDomain.show()

#construction of the Mesh--------------------------------------
MyMesh = Mesh()
MyMesh.generate_mesh_from_domain(MyDomain, 0.02)
MyMesh.show()
MyMesh.save_to_json(str(Path(__file__).parent / "data" / "RuGrandmont"))


#Construction of the initial datum ---------------------------------------
MyMap = CellValueMap(MyMesh)
MyMap.generate_random()
MyMap.show()


#Construction of the solver and computing until the domain is empty --------------
opt = dict(model = "hughes",
            filename = str(Path(__file__).parent / "data" / "RuGrandmont"),
            save = True,
            verbose = True, #Set to false for less information printed in the console
            lwrSolver = {   'convexFlux' : True,
                            'anNum' : "dichotomy",
                            'method' : "midVector",
                            'ApproximationThreshold' : 0.0001},
            eikoSolver = {  'constrained' : True,
                            'NarrowBandDepth' : 2})

Solver = PedestrianSolver(MyMesh, 0.01, initial_density = MyMap, options=opt)
Solver.compute_until_empty(5000)
