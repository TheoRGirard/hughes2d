
import json
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.animation as animation

import csv

filename = "../../../Visualiseur/data/examples/additional_computations_test3"

values = []
numFrames = 0

with open(filename +"_zones_mean_densiy.csv", mode ='r')as file:
    csvFile = csv.reader(file)
    numLine = 0
    for lines in csvFile:
        #print("ligne de longueur : ", len(lines))
        if(numLine == 0):
            values = [[] for i in range(1,len(lines))]
            numFrames = len(lines)-1
        else:
            for i in range(len(lines)-1):
                #print(i,"/", len(values))
                values[i].append(float(lines[i+1]))
        numLine += 1

#print(values)

fig, ax = plt.subplots()
xdata, ydata = np.arange(0,5,0.025), np.zeros(len(np.arange(0,5,0.025)))
print(len(values[1]),"/",len(np.arange(0,5,0.1)))
ln = ax.fill_between(xdata, ydata, ydata, color="darkslategrey")

def init():
    ax.set_xlim(0, 5)
    ax.set_ylim(0.0, 1.0)
    ln = ax.fill_between(xdata, ydata, ydata, color="darkslategrey")
    return ln,

def update(frame):
    ln = ax.fill_between(xdata,ydata,values[frame], color="darkslategrey")
    return ln,

ani = animation.FuncAnimation(fig, update, frames=range(numFrames), init_func=init, blit=True, interval=100,repeat=True)
plt.show()
