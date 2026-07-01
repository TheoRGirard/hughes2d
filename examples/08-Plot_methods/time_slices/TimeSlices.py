from pathlib import Path

from hughes2d import Plotter

Plotter.save_time_slices(
    [0, 1, 2],
    filename=str(Path(__file__).parent / "time_slices"),
    slicename=str(Path(__file__).parent.parent / "fig" / "slice"),
    limits=[[-2, 12], [0, 7]],
    ft_size=10,
    dpi_set=1000,
)
