"""
© Copyright 2025 Girard Théo

This file is part of the Hughes2d package.

The Hughes2d package is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.

The Hughes2d package is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the
Hughes2d package. If not, see <https://www.gnu.org/licenses/>.
"""

import csv
import json
from pathlib import Path

import numpy as np

FigureType = None

try:
    import matplotlib.pyplot as plt
    from matplotlib import animation, cm, collections
    #FigureType = FigureType | plt.Figure
except ImportError:
    plt = None

try:
    import plotly.figure_factory as ff
    import plotly.graph_objects as go
    #FigureType = FigureType | go.Figure
except ImportError:
    go, ff = None, None



EMPTY_LIST = []

def convert_to_mp4(filename:str, limits:list[list[float]] = EMPTY_LIST, dpi_set:int = 300, pic_size:tuple[int,int] = (8,6), *, plot_scale:bool = True) -> None:
    """Create a video file from the data files of a simulation.

    Args:
        filename (str): the path and basename for the data files.
        limits (list[list[float]], optional): the scope of the simulation as
            [[x_min,x_max],[y_min,y_max]].
        dpi_set (int, optional): resolution of the video as dots per inches (dpi)
        pic_size (tuple[int,int], optional): size of the picture as (width,height)
        plot_scale (bool, keyword-only): toggles the plot of a color scale for
            the density.

    Raises:
        ImportError: if matplotlib is not installed.

    """
    if not plt:
        msg = "matplotlib is required for this function."
        raise ImportError(msg)

    with Path(filename +"_mesh.json").open("r") as f:
        data = json.load(f)
        triangles = np.array([ [ [ data["vertices"][i][0], data["vertices"][i][1] ]
                                for i in triangle] for triangle in data["triangles"] ])


    values = []

    with Path(filename +"_densities.csv").open("r") as file:
        csv_file = csv.reader(file)
        for lines in csv_file:
            values = [*values,
                      np.array([float(lines[i]) for i in range(len(lines))])]


    fig, ax = plt.subplots(figsize= pic_size, dpi=dpi_set)

    col = collections.PolyCollection(triangles)
    col.set_cmap(cm.viridis)
    col.set_clim([0, 1])
    rgcol = col.set_array(values[0])

    ax.add_collection(col)
    if plot_scale:
        fig.colorbar(rgcol, ax=ax, label="density")


    if len(limits) > 0:
        ax.set_xlim(limits[0][0],limits[0][1])
        ax.set_ylim(limits[1][0],limits[1][1])
        plt.axis("equal")


    def update(frame:int):
        col.set_array(values[frame])
        ax.add_collection(col)

        return (col,ax)


    ani = animation.FuncAnimation(fig=fig, func=update, frames=len(values), interval=30)

    ff_writer = animation.FFMpegWriter(fps=25)
    ani.save(filename+".mp4", writer = ff_writer, dpi=dpi_set)


def save_time_slices(times:list[float], filename:str, slicename:str, limits:list[list[float]] = EMPTY_LIST, dpi_set:int = 300, pic_size: tuple[int,int] = (8,6), *, plot_scale:bool = True) -> None:
    """Export an image file from the data files of a simulation for each time slice required in the ``times`` parameter.

    Args:
        times (list[float]): a list of the timing (in seconds) when a picture should be
            extracted from the simulation.
        filename (str): the path and basename for the data files.
        slicename (str): the path and basename for the exported pictures.
        limits (list[list[float]]): the scope of the simulation as
            [[x_min,x_max],[y_min,y_max]].
        dpi_set (int): resolution of the video as dots per inches (dpi)
        pic_size (tuple[int,int], optional): size of the picture as (width,height)
        plot_scale (bool, keyword-only): toggles the plot of a color scale for
            the density.

    Raises:
        ImportError: if matplotlib is not installed.

    """
    if not plt:
        msg = "matplotlib is required for this function."
        raise ImportError(msg)

    int_times = [int(times[i]*25) for i in range(len(times))]

    with Path(filename +"_mesh.json").open("r") as f:
        data = json.load(f)
        triangles = np.array([ [ [ data["vertices"][i][0], data["vertices"][i][1] ]
                                for i in triangle] for triangle in data["triangles"] ])


    values = []

    with Path(filename +"_densities.csv").open("r") as file:
        csv_file = csv.reader(file)
        for lines in csv_file:
            values = [*values,
                      np.array([float(lines[i]) for i in range(len(lines))])]

    for i,t in enumerate(times):

        if int_times[i] < len(values):
            fig, ax = plt.subplots(figsize= pic_size, dpi=dpi_set)

            col = collections.PolyCollection(triangles)
            col.set_cmap(cm.viridis)
            col.set_clim([0, 1])
            rgcol = col.set_array(values[int_times[i]])


            if len(limits) > 0:
                ax.set_xlim(limits[0][0],limits[0][1])
                ax.set_ylim(limits[1][0],limits[1][1])
                plt.axis("equal")

            ax.add_collection(col)
            ax.set_title("t = "+str(t)+"s")
            if plot_scale:
                fig.colorbar(rgcol, ax=ax, label="density")
            plt.savefig(slicename + str(t) +"s.png", dpi=dpi_set)
        else:
            print("Warning: time slice ", t,
                  " ignored because the simulation is too short.")



def plot_vector_field(vertex_list:list[list[float]], triangle_list:list[list[int]], vector_field:list[list[float]], *, plot_mesh:bool=True) -> None:
    """Display a vector field passed as a parameter with plotly.

    Args:
        vertex_list (list[list[float]]): a list of the vertices of the mesh on which
            the vector field will be plotted.
        triangle_list (list[list[int]]): a list of the triangles of the mesh on which
            the vector field will be plotted.
        vector_field (list[list[float]]): the vector field to plot. Must be of the same
            shape as ``triangle_list``.
        plot_mesh (bool, keyword-only): whether or not the mesh should be plotted in
            background.

    Raises:
        ImportError: if plotly is not installed.

    """
    if not go:
        msg = "plotly is required for this function."
        raise ImportError(msg)

    fig = go.Figure()
    if(plot_mesh):
        for triangle in triangle_list:
            fig.add_trace(go.Scatter(x=([vertex_list[i][0] for i in triangle]
                                        +[vertex_list[triangle[0]][0]]),
                                    y=([vertex_list[i][1] for i in triangle]
                                       +[vertex_list[triangle[0]][1]]),
                            fill="toself",
                            fillcolor="White",
                            mode="lines",
                            line={
                                "color" : "Black",
                                "width" : 1,
                             }))
        fig.update_layout(yaxis={
                "scaleanchor" : "x",
                "scaleratio" : 1,
                })
    fig_quiv = ff.create_quiver([(vertex_list[triangle[0]][0]+vertex_list[triangle[1]][0]
                                 +vertex_list[triangle[2]][0])/3 for triangle in triangle_list],
                                [(vertex_list[triangle[0]][1]+vertex_list[triangle[1]][1]
                                  +vertex_list[triangle[2]][1])/3 for triangle in triangle_list],
                                [V[0] for V in vector_field], [V[1] for V in vector_field])
    fig.add_traces(fig_quiv.data)
    fig.update_layout(yaxis={
            "scaleanchor" : "x",
            "scaleratio" : 1,
            })
    fig.show()

def add_vector_field_plot(fig:FigureType, vertex_list:list[float], triangle_list:list[list[int]], vector_field:list[list[float]], color:list[int]=[255,255,255,0.9]) -> None:
    """Add a vector field plot to a given figure.

    Args:
        fig (go.Figure): the figure on which the vector field plot should be added;
        vertex_list (list[list[float]]): a list of the vertices of the mesh on which
            the vector field will be plotted.
        triangle_list (list[list[int]]): a list of the triangles of the mesh on which
            the vector field will be plotted.
        vector_field (list[list[float]]): the vector field to plot. Must be of the same
            shape as ``triangle_list``.
        plot_mesh (bool, keyword-only): whether or not the mesh should be plotted in
            background.
        color (list[int], optional): the color of the vectors.

    Raises:
        ImportError: if plotly is not installed.

    """
    if not go:
        msg = "plotly is required for this function."
        raise ImportError(msg)

    fig_quiv = ff.create_quiver([(vertex_list[triangle[0]][0]+vertex_list[triangle[1]][0]
                                  +vertex_list[triangle[2]][0])/3
                                 for triangle in triangle_list],
                                [(vertex_list[triangle[0]][1]+vertex_list[triangle[1]][1]
                                  +vertex_list[triangle[2]][1])/3
                                 for triangle in triangle_list],
                                [V[0] for V in vector_field],
                                [V[1] for V in vector_field],
                                line={
                                    "color":"rgba("+str(color[0])+","+str(color[1])+","+str(color[2])+","+str(color[3])+")",
                                })
    fig.add_traces(fig_quiv.data)
    fig.update_layout(yaxis={
            "scaleanchor" : "x",
            "scaleratio" : 1,
            })
