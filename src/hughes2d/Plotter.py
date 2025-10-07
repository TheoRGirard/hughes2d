
import json
import numpy as np
import csv

try:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.collections as collections
    import matplotlib.animation as animation
except ImportError:
    plt = None

try:
    import plotly.figure_factory as ff
    import plotly.graph_objects as go
except ImportError:
    go, ff = None, None

def convertToMP4(filename:str, limits:List[List[float]] = [], dpi_set:int = 300) -> None:
    """
    Creates a video file from the data files of a simulation.

    Args:
        filename (str): the path and basename for the data files.
        limits (List[List[float]]): the scope of the simulation as [[x_min,x_max],[y_min,y_max]].
        dpi_set (int): resolution of the video as dots per inches (dpi)

    Raises:
        ImportError: if matplotlib is not installed.
    """

    if not matplotlib:
        raise ImportError("matplotlib is required for this function.")

    with open(filename +"_mesh.json") as f:
        data = json.load(f)
        Triangles = np.array([ [ [ data['vertices'][i][0], data['vertices'][i][1] ] for i in triangle] for triangle in data['triangles'] ])


    values = []

    with open(filename +"_densities.csv", mode ='r')as file:
        csvFile = csv.reader(file)
        for lines in csvFile:
            values.append(np.array([float(lines[i]) for i in range(len(lines))]))


    fig, ax = plt.subplots(dpi=dpi_set)

    col = collections.PolyCollection(Triangles)
    col.set_cmap(cm.viridis)
    col.set_clim([0, 1])
    rgcol = col.set_array(values[0])

    ax.add_collection(col)
    fig.colorbar(rgcol, ax=ax, label="density")


    if len(limits) > 0:
        ax.set_xlim(limits[0][0],limits[0][1])
        ax.set_ylim(limits[1][0],limits[1][1])
        plt.axis('equal')


    def update(frame):
        # for each frame, update the data stored on each artist.
        col.set_array(values[frame])
        ax.add_collection(col)

        return (col,ax)


    ani = animation.FuncAnimation(fig=fig, func=update, frames=len(values), interval=30)

    FFwriter = animation.FFMpegWriter(fps=25)
    ani.save(filename+'.mp4', writer = FFwriter, dpi=dpi_set)


def saveTimeSlices(times:List[float], filename:str, slicename:str, limits:List[List[float]] = [], dpi_set:int = 300) -> None:
    """
    Exports an image file from the data files of a simulation for each time slice required in the ``times`` parameter.

    Args:
        times (List[float]): a list of the timing (in seconds) when a picture should be extracted from the simulation.
        filename (str): the path and basename for the data files.
        slicename (str): the path and basename for the exported pictures.
        limits (List[List[float]]): the scope of the simulation as [[x_min,x_max],[y_min,y_max]].
        dpi_set (int): resolution of the video as dots per inches (dpi)

    Raises:
        ImportError: if matplotlib is not installed.
    """
    if not matplotlib:
        raise ImportError("matplotlib is required for this function.")

    IntTimes = [int(times[i]*25) for i in range(len(times))]

    with open(filename +"_mesh.json") as f:
        data = json.load(f)
        Triangles = np.array([ [ [ data['vertices'][i][0], data['vertices'][i][1] ] for i in triangle] for triangle in data['triangles'] ])


    values = []

    with open(filename +"_densities.csv", mode ='r')as file:
        csvFile = csv.reader(file)
        for lines in csvFile:
            values.append(np.array([float(lines[i]) for i in range(len(lines))]))


    for i,t in enumerate(times):
        fig, ax = plt.subplots()

        col = collections.PolyCollection(Triangles)
        col.set_cmap(cm.viridis)
        col.set_clim([0, 1])
        rgcol = col.set_array(values[IntTimes[i]])


        if len(limits) > 0:
            ax.set_xlim(limits[0][0],limits[0][1])
            ax.set_ylim(limits[1][0],limits[1][1])
            plt.axis('equal')

        ax.add_collection(col)
        ax.set_title("t = "+str(times[i])+"s")
        fig.colorbar(rgcol, ax=ax, label="density")
        plt.savefig(slicename + str(times[i]) +"s.png", dpi=dpi_set)



def plotVectorField(VertexList:List[List[float]], TriangleList:List[List[int]], VectorField:List[List[float]], plotMesh:bool=True) -> None:
    """
    Displays a vector field passed as a parameter with plotly.

    Args:
        VertexList (List[List[float]]): a list of the vertices of the mesh on which the vector field will be plotted.
        TriangleList (List[List[int]]): a list of the triangles of the mesh on which the vector field will be plotted.
        VectorField (List[List[float]]): the vector field to plot. Must be of the same shape as ``TriangleList``.
        plotMesh (bool, optional): whether or not the mesh should be plotted in background.

    Raises:
        ImportError: if plotly is not installed.
    """
    if not go:
        raise ImportError("plotly is required for this function.")

    fig = go.Figure()
    if(plotMesh):
        for T in TriangleList:
            fig.add_trace(go.Scatter(x=[VertexList[i][0] for i in T]+[VertexList[T[0]][0]],
                                    y=[VertexList[i][1] for i in T]+[VertexList[T[0]][1]],
                            fill="toself",
                            fillcolor="White",
                            mode="lines",
                            line=dict(
                                color="Black",
                                width=1
                             )))
        fig.update_layout(yaxis=dict(
            scaleanchor='x',
            scaleratio=1))
    figQuiv = ff.create_quiver([(VertexList[T[0]][0]+VertexList[T[1]][0]+VertexList[T[2]][0])/3 for T in TriangleList],
                                [(VertexList[T[0]][1]+VertexList[T[1]][1]+VertexList[T[2]][1])/3 for T in TriangleList],
                                [V[0] for V in VectorField], [V[1] for V in VectorField])
    fig.add_traces(figQuiv.data)
    fig.update_layout(yaxis=dict(
        scaleanchor='x',
        scaleratio=1))
    fig.show()

def addVectorFieldPlot(fig, VertexList, TriangleList, VectorField, color=[255,255,255,0.9]):
    if not go:
        raise ImportError("plotly is required for this function.")

    figQuiv = ff.create_quiver([(VertexList[T[0]][0]+VertexList[T[1]][0]+VertexList[T[2]][0])/3 for T in TriangleList],
                                [(VertexList[T[0]][1]+VertexList[T[1]][1]+VertexList[T[2]][1])/3 for T in TriangleList],
                                [V[0] for V in VectorField], [V[1] for V in VectorField],
                                line=dict(
                                    color='rgba('+str(color[0])+','+str(color[1])+','+str(color[2])+','+str(color[3])+')')
                                )
    fig.add_traces(figQuiv.data)
    fig.update_layout(yaxis=dict(
        scaleanchor='x',
        scaleratio=1))
