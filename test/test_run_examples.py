import runpy
from pathlib import Path

import pytest

try:
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

try:
    import plotly.graph_objects as go
except ImportError:
    go = None

example_folder = Path(__file__).parent.parent / "examples"

example_files_with_plot = [Path(example_folder, folder, file)
                 for folder, file in [
          ("02-Simulating_Hughes_1", "Simulating_Grandmont_Restaurant.py"),
          ("03-Simulating_Hughes_2", "Comparison_Hughes_vs_explicit_vector_field.py"),
          ("04-Simulating_Eikonal_equation", "Simulating_shortest_path.py"),
          ("05-Simulating_LWR", "Simulating_LWR.py"),
                                   ]]

example_files_without_plot = [Path(example_folder, folder, file)
                 for folder, file in [
          ("00-getting_started","getting_started.py"),
          ("01-Mesh_generation", "Generating_square_mesh.py"),
          ("07-Additional_Computations", "AdditionalComputations.py"),
                                   ]]

@pytest.mark.parametrize("filename", example_files_without_plot)
def test_run_basic_examples(filename):
    runpy.run_path(filename)

@pytest.mark.parametrize("filename", example_files_with_plot)
def test_run_plot_examples(filename):
    if plt or go:
        runpy.run_path(filename)
    else:
        with pytest.raises(ImportError):
            runpy.run_path(filename)
