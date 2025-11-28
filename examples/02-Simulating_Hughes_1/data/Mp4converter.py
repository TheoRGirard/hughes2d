from pathlib import Path

from hughes2d import Plotter

filename = str(Path(__file__).parent / "RuGrandmont")
Plotter.convertToMP4(filename, limits=[[1,13],[-11,-3]])
