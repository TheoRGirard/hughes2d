
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.collections as collections

import matplotlib.animation as animation

import plotly.figure_factory as ff
import plotly.graph_objects as go

import csv

def plot(filename):

    with open(filename +"_mesh.json") as f:
        data = json.load(f)
        Triangles = np.array([ [ [ data['vertices'][i][0], data['vertices'][i][1] ] for i in triangle] for triangle in data['triangles'] ])


    values = []

    with open(filename +"_densities.csv", mode ='r')as file:
        csvFile = csv.reader(file)
        for lines in csvFile:
            values.append(np.array([float(lines[i]) for i in range(len(lines))]))


    fig, ax = plt.subplots()

    col = collections.PolyCollection(Triangles)
    col.set_array(values[0])
    col.set_cmap(cm.viridis)
    ax.add_collection(col)

    ax.set_xlim(0,12)
    ax.set_ylim(0,7)

    def update(frame):
        # for each frame, update the data stored on each artist.
        col.set_array(values[frame])
        col.set_cmap(cm.viridis)
        ax.add_collection(col)
        return (col,ax)


    ani = animation.FuncAnimation(fig=fig, func=update, frames=len(values), interval=30)

    plt.show()

def plotVectorField(VertexList, TriangleList, VectorField):
    if(go):
        fig = go.Figure()
        #self.Mesh.domain.addPlot(fig)
        figQuiv = ff.create_quiver([(VertexList[T[0]][0]+VertexList[T[1]][0]+VertexList[T[2]][0])/3 for T in TriangleList],
                                    [(VertexList[T[0]][1]+VertexList[T[1]][1]+VertexList[T[2]][1])/3 for T in TriangleList],
                                    [V[0] for V in VectorField], [V[1] for V in VectorField])
        fig.add_traces(figQuiv.data)
        fig.update_layout(yaxis=dict(
            scaleanchor='x',
            scaleratio=1))
        fig.show()
    else:
        raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")
