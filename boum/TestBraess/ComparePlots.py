import csv
import matplotlib.pyplot as plt

import numpy as np

directory = "data/"
files = ["Hughes"]+["HughesX" + str(np.round(x,2)) for x in np.arange(3,4.1,0.1) ] #np.arange(3,3.09,0.1)
plotnames = ["_total_mass","_max_density","_zones_mean_densiy"]

modeZone = "multCurves"
numZones = 2

Threshold = 0.1
TimeForThreshold = []

valueAt15 = []
valueAt20 = []
Slopes = []

foundT,found15,found20 = False,False,False

for plotName in plotnames:
    if(plotName == "_zones_mean_densiy"):
        if(modeZone == "multCurves"):
            fig, axs = plt.subplots(numZones)
            fig.suptitle('Mean density in different zones')
            for model in files:
                with open(directory+model+plotName+".csv", mode ='r')as file:
                    csvFile = csv.reader(file)
                    numLine = 0

                    zoneData = []
                    zoneNames = []
                    for lines in csvFile:
                        if(numLine == 0):
                            Times = [float(lines[i]) for i in range(1,len(lines))]
                        else:
                            zoneData.append([float(lines[i]) for i in range(1,len(lines))])
                            zoneNames.append(lines[0])
                        numLine += 1
                    for i in range(numZones):
                        axs[i].plot(Times,zoneData[i], label = model)
                        axs[i].set_title(zoneNames[i])
                        axs[i].legend()
            plt.show()

    else:
        fig, ax = plt.subplots()

        for model in files:
            foundT,found15,found20 = False,False,False
            with open(directory+model+plotName+".csv", mode ='r')as file:
                csvFile = csv.reader(file)
                numLine = 0
                for lines in csvFile:
                    if(numLine == 0):
                        Times = [float(lines[i]) for i in range(len(lines))]
                    else:
                        Data = [float(lines[i]) for i in range(len(lines))]
                        if(plotName == "_total_mass"):
                            for i in range(len(lines)):
                                if(Data[i] <= Threshold) and not foundT:
                                    TimeForThreshold.append(Times[i])
                                    foundT = True
                                if(Times[i] >= 15) and not found15:
                                    valueAt15.append(Data[i])
                                    found15 = True
                                if(Times[i] >= 20) and not found20:
                                    valueAt20.append(Data[i])
                                    found20 = True
                    numLine += 1

                ax.plot(Times, Data, label = model)
                ax.set_title(plotName + " with respect to time")
        ax.legend()
        plt.show()

print(TimeForThreshold)

fig, ax = plt.subplots()
ax.plot(np.arange(3,4.4,0.1), [TimeForThreshold[0] for _ in np.arange(3,4.4,0.1)],"r--", label = "evac. time without an obstacle")
ax.plot(np.arange(3,4.1,0.1), TimeForThreshold[1:], label = "evacuation time with obstacle")
ax.set_title("Evacuation time with respect to the position of the obstacle.")
ax.set_xlabel("Xc")
ax.set_ylabel("Evac. time (s)")
ax.legend()
plt.show()

for i in range(len(valueAt15)):
    Slopes.append((valueAt15[i]-valueAt20[i])/5)

print(Slopes)
print(valueAt15)
print(valueAt20)

fig, ax = plt.subplots()
ax.plot(np.arange(3,4.4,0.1), [Slopes[0] for _ in np.arange(3,4.4,0.1)],"r--", label = "evacuation speed without pillar")
ax.plot(np.arange(3,4.1,0.1), Slopes[1:], label = "time of evacuation at 0.1")
ax.set_title("Evacuation speed with respect to the position of the pillar")
ax.legend()
plt.show()
