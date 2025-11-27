import pytest

from hughes2d.Mesh2D import Mesh, NonConvexDomain

"""
This file tests the import and export methods for the *Mesh* and *NonConvexDomain*
classes.
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
        Domain1.import_from_dxf("test/ressources/test_Domain.dxf")
        Domain1.show()
        assert [200,100] in Domain1
        assert [-3,-3] not in Domain1
    else:
        with pytest.raises(ImportError):
            Domain1.import_from_dxf("ressources/test_Domain.dxf")

MyMesh = Mesh()

def test_import_msh():
    if meshio:
        MyMesh.import_mesh_from_msh("test/ressources/test_mesh.msh")
        assert len(MyMesh.triangles) == 111
    else:
        with pytest.raises(ImportError):
            MyMesh.import_mesh_from_msh("filename.msh")

def test_import_msh_FF():
    if meshio:
        MyMesh.import_mesh_from_msh_free_fem("test/ressources/test_mesh_FreeFEM.msh", flag_dict = {"domain" : 0, "exit" : [2], "wall" : [1,3,4]})
        assert len(MyMesh.triangles) == 72
    else:
        with pytest.raises(ImportError):
            MyMesh.import_mesh_from_msh("filename.msh")
