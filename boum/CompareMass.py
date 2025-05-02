import csv
import matplotlib.pyplot as plt

file1 = "LargeurCouloir0.2_total_mass.csv"
file2 = "LargeurCouloir0.8_total_mass.csv"

with open(file1, mode ='r')as file:
    csvFile = csv.reader(file)
    numLine = 0
    for lines in csvFile:
        if(numLine == 0):
            Times = [float(lines[i]) for i in range(len(lines))]
        elif(numLine == 1):
            Mass1 = [float(lines[i]) for i in range(len(lines))]
        numLine += 1

with open(file2, mode ='r')as file:
    csvFile = csv.reader(file)
    numLine = 0
    for lines in csvFile:
        if(numLine == 0):
            Times2 = [float(lines[i]) for i in range(len(lines))]
        elif(numLine == 1):
            Mass2 = [float(lines[i]) for i in range(len(lines))]
        numLine += 1

plt.plot(Times, Mass1)
plt.plot(Times2, Mass2)
plt.show()
