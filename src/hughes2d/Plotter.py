
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.collections as collections

import matplotlib.animation as animation

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
