import runpy
from pathlib import Path

import pytest

try:
    import matplotlib as mpl
    mpl.use("Agg")
except ImportError:
    pass

example_folder = Path(__file__).parent.parent / "examples"

example_files = [Path(example_folder, folder, file)
                 for folder, file in [
          ("00-getting_started","getting_started.py"),
          ("01-Mesh_generation", "Generating_square_mesh.py"),
          ("02-Simulating_Hughes_1", "Simulating_Grandmont_Restaurant.py"),
          ("03-Simulating_Hughes_2", "Comparison_Hughes_vs_explicit_vector_field.py"),
          ("04-Simulating_Eikonal_equation", "Simulating_shortest_path.py"),
          ("05-Simulating_LWR", "Simulating_LWR.py"),
          ("07-Additional_Computations", "AdditionalComputations.py"),
                                   ]]

@pytest.mark.parametrize("filename", example_files)
def test_run_basic_examples(filename):
    runpy.run_path(filename)

def test_run_second_exec():
    runpy.run_path(Path(__file__).parent.parent / "examples"
                   / "01-Mesh_generation"/ "Generating_square_mesh.py")
