import numpy as np
from hughes2d.Plotter import *



ListModel = ["data/TestTPID0"] + [ "data/TestTPID0_"+str(i) for i in range(3,7)]

TM = np.arange(16,0,-0.1)
Indices = [[] for _ in ListModel]
Masses = [[] for _ in ListModel]
Seuil = 0.02
T = [[] for _ in ListModel]
Data = []
Times = []

for i in range(len(ListModel)):
    with open(ListModel[i] +"_total_mass.csv", mode ='r') as file:
        csvFile = csv.reader(file)
        numLine = 0
        for lines in csvFile:
            if(numLine == 0):
                Times = [float(lines[i]) for i in range(len(lines))]
            else:
                Data = [float(lines[i]) for i in range(len(lines))]
            numLine += 1
        print(len(Data))
        j = 0
        Numass = 0
        while Numass < len(TM):
            if(j >= len(Data)):
                break
            if(TM[Numass] > Data[j] + Seuil):
                if(j == 0):
                    print("Adding NaN")
                    Indices[i].append("NaN")
                else:
                    Indices[i].append(max(0,j-1))
                Masses[i].append(TM[Numass])
                Numass += 1
            else:
                if abs(TM[Numass] - Data[j]) < Seuil:
                    Indices[i].append(j)
                    Masses[i].append(TM[Numass])
                    if(abs(TM[Numass] - 13.5) < Seuil):
                        T[i].append(Times[j])
                    if(abs(TM[Numass] - 14.5) < Seuil):
                        T[i].append(Times[j])
                    Numass += 1
                else:
                    j += 1


if len(Masses[0]) != len(Masses[1]) or len(Masses[0]) != len(Masses[2]) :
    print([len(Masses[i]) for i in range(3)])
    print(len(TM))
    raise ValueError("ça marche po")

i = 0
#saveTimeSlices([0] + T[i],filename=ListModel[i],slicename = "figs/testTP_ID0_TM", limits=[[0,12],[0,7]])



Datas = [[] for _ in ListModel]
TrianglesAreas = []
with open("testTP" +"_mesh.json") as f:
    data = json.load(f)
    TrianglesAreas = np.array(data['cellAreas'])



for i in range(len(ListModel)):
    with open(ListModel[i]+"_densities.csv", mode ='r')as file:
        csvFile = csv.reader(file)
        numLine = 0
        for lines in csvFile:
            if numLine in Indices[i]:
                Datas[i].append(np.array([float(lines[j]) for j in range(len(lines))]))
            numLine += 1

for i in range(len(ListModel)):
    k = 0
    while( Indices[i][k] == "NaN"):
        k += 1
    print(k)
    Datas[i] = ["Nan" for _ in range(k)] + Datas[i]
print([len(Datas[i]) for i in range(3)])

L1Diffs = []
X = []


for j,m in enumerate(TM):
    max = 0
    for num1 in range(len(ListModel)):
        for num2 in range(num1 + 1, len(ListModel)):
            if(Indices[num1][j] != "NaN" and Indices[num2][j] != "NaN"):
                D = 0
                for k in range(len(TrianglesAreas)):
                    D += abs(Datas[num1][j][k]-Datas[num2][j][k])*TrianglesAreas[k]
                if(D > max):
                    max = D

    L1Diffs.append(max/72)
    X.append(m)




fig, axs = plt.subplots(1)

axs.plot(X,L1Diffs)
axs.set_xlabel("total mass")
axs.set_xlim(16,0)
axs.set_ylabel("Diff L¹")
#plt.show()
plt.savefig("figs/comparaisonID0-random.png")
