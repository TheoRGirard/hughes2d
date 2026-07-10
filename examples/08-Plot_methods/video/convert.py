from pathlib import Path

import hughes2d.Plotter

filename = str(Path(__file__).parent / "Hughes")

hughes2d.Plotter.convert_to_mp4(filename)
