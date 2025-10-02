from hughes2d.Mesh2D import *
import pytest

"""
This file tests the import and export methods for the *Mesh* and *NonConvexDomain* classes.
"""

try:
    import ezdxf
except ImportError:
    ezdxf = None

try:
    import meshio
except ImportError:
    meshio = None

Domain1 = NonConvexDomain()

def test_import_dxf_domain():
    if ezdxf:
        Domain1.importFromDXF("test/ressources/test_Domain.dxf")
        Domain1.show()
        assert [200,100] in Domain1
        assert [-3,-3] not in Domain1
    else:
        with pytest.raises(ImportError):
            Domain1.importFromDXF("ressources/test_Domain.dxf")

MyMesh = Mesh()

def test_import_msh():
    if meshio:
        MyMesh.importMeshFromMsh("test/ressources/test_mesh.msh")
        assert len(MyMesh.triangles) == 111
    else:
        with pytest.raises(ImportError):
            MyMesh.importMeshFromMsh("filename.msh")

def test_import_msh_FF():
    if meshio:
        MyMesh.importMeshFromMshFreeFem("test/ressources/test_mesh_FreeFEM.msh")
        assert len(MyMesh.triangles) == 5038
    else:
        with pytest.raises(ImportError):
            MyMesh.importMeshFromMsh("filename.msh")
