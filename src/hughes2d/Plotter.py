
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.collections as collections

import matplotlib.animation as animation

import plotly.figure_factory as ff
import plotly.graph_objects as go

import csv

def convertToMP4(filename):

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
    #fig.colorbar(ax,label="density")
    ax.set_xlim(0,12)
    ax.set_ylim(0,7)


    def update(frame):
        # for each frame, update the data stored on each artist.
        col.set_array(values[frame])
        col.set_cmap(cm.viridis)
        ax.add_collection(col)
        return (col,ax)


    ani = animation.FuncAnimation(fig=fig, func=update, frames=len(values), interval=30)

    FFwriter = animation.FFMpegWriter(fps=25)
    ani.save(filename+'.mp4', writer = FFwriter)


def saveTimeSlices(times, filename, slicename, limits = []):

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
        rgcol = col.set_array(values[IntTimes[i]])
        col.set_cmap(cm.viridis)

        if len(limits) > 0:
            ax.set_xlim(limits[0][0],limits[0][1])
            ax.set_ylim(limits[1][0],limits[1][1])
            plt.axis('equal')

        ax.add_collection(col)
        ax.set_title("t = "+str(times[i])+"s")
        fig.colorbar(rgcol, ax=ax, label="density")
        plt.savefig("fig/" + slicename + str(times[i]) +"s.png")



def plotVectorField(VertexList, TriangleList, VectorField, plotMesh=True):
    if(go):
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
    else:
        raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")

def addVectorFieldPlot(fig, VertexList, TriangleList, VectorField, color=[255,255,255,0.9]):
    if(go):
        #self.Mesh.domain.addPlot(fig)
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
    else:
        raise ImportError("No plotting module found. Try installing plotly or matplotlib if you want to use show methods")
