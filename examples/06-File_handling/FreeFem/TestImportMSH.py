from pathlib import Path

from hughes2d import Mesh

filename = str(Path(__file__).parent / "mesh_FF.msh")

MyMesh = Mesh()
MyMesh.import_mesh_from_msh_free_fem(filename)
MyMesh.show()
