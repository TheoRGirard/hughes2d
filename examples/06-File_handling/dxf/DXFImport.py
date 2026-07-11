from pathlib import Path

from hughes2d.Mesh2D import Mesh, NonConvexDomain

filename = Path(__file__).parent / "config_simple.dxf"

MyDomain = NonConvexDomain()
MyDomain.import_from_dxf(filename)
MyDomain.show()

MyMesh = Mesh()
MyMesh.generate_mesh_from_domain(MyDomain, 10)
MyMesh.show()
