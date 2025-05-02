import csv
import matplotlib.pyplot as plt

directory = "data/"
#files = ["HughesID1","LWRconstantFieldID1","ColomboGaravelloID1"]
files = ["HughesID2","LWRconstantFieldID2","ColomboGaravelloID2"]
plotnames = ["_total_mass","_max_density","_zones_mean_densiy"]

modeZone = "multCurves"
numZones = 4
zoneNames2 = [ "right side of the room", "left corridor","left side of the room",  "right corridor"]

for plotName in plotnames:
    if(plotName == "_zones_mean_densiy"):
        if(modeZone == "multCurves"):
            #fig, axs = plt.subplots(numZones)
            #fig.suptitle('Mean density in different zones')
            for i in range(numZones):
                fig, axs = plt.subplots(1)
                fig.suptitle('Mean density in the '+ zoneNames2[i])
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

                        axs.plot(Times,zoneData[i], label = model[:-3])
                        #axs[i].set_title(zoneNames[i])
                        axs.set_xlabel("t")
                        axs.legend()
                plt.savefig("figs/meanDensity"+zoneNames[i]+".png")
            #plt.show()

    else:
        fig, ax = plt.subplots()
        for model in files:
            with open(directory+model+plotName+".csv", mode ='r')as file:
                csvFile = csv.reader(file)
                numLine = 0
                for lines in csvFile:
                    if(numLine == 0):
                        Times = [float(lines[i]) for i in range(len(lines))]
                    else:
                        Data = [float(lines[i]) for i in range(len(lines))]
                    numLine += 1

                ax.plot(Times, Data, label = model[:-3])
                ax.set_xlabel("t")
                ax.set_ylabel("tot. mass")
                ax.legend()
                ax.set_title("Total mass " + "with respect to time")
        ax.legend()
        plt.show()
