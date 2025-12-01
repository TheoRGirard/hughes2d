from pathlib import Path

import hughes2d.Plotter

filename = str(Path(__file__).parent / "Hughes")

hughes2d.Plotter.convertToMP4(filename)
