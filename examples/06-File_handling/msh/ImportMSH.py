from pathlib import Path

from hughes2d import Mesh

#From : https://github.com/nschloe/meshio/blob/main/tests/meshes/msh/insulated-2.2.msh
filename = str(Path(__file__).parent / "insulated-2.2.msh")


MyMesh = Mesh()
MyMesh.import_mesh_from_msh(filename)
MyMesh.show()
