from hughes2d import *
import random as alea
import numpy as np
import matplotlib.pyplot as plt

MyDomain = Mesh2D.NonConvexDomain()
MyDomain = Mesh2D.NonConvexDomain([[-1,0],[-1,2],[1,2],[1,0]])
MyDomain.addExit([[-1,0],[1,0]])

MaxDiffs = [[],[],[]] #order TVE vs explicit, TVTuD2 vs explicit, TVTC vs explicit
Means = [[],[],[]]
VectorDiffMax = [[],[],[]]
VectorDiffMeans = [[],[],[]]
VectorDiffMaxE = [[],[],[]]
VectorDiffMeansE = [[],[],[]]

X = np.arange(100, 1100, 100)

for d in X:
    print("Step : ", 1/d)
    MyMesh = Mesh2D.Mesh()
    MyMesh.generateMeshFromDomain(MyDomain, 1/d)
    #MyMesh.show()

    MyMap = Mesh2D.CellValueMap(MyMesh)

    for index, triangle in enumerate(MyMesh.triangles):
        MyMap[index] = 1
        for pair in [[0,1],[1,2],[2,0]]:
            if(MyMesh.vertices[triangle[pair[0]]][0])*(MyMesh.vertices[triangle[pair[1]]][0]) <= 0:
                MyMap.values[index] = 0

    ConstantVectorField = []

    for P in MyMesh.barycenters:
        if(3*abs(P[0])/np.sqrt(3)>=P[1]):
            ConstantVectorField.append([0,-1])
        elif(abs(P[0]) < 0.001):
            ConstantVectorField.append([0,-1])
        else:
            V = [-P[0],abs(P[0]/3)]
            ConstantVectorField.append([V[0]/np.sqrt(V[0]**2 + V[1]**2), V[1]/np.sqrt(V[0]**2 + V[1]**2)])

    ExplicitSolution = []

    for P in MyMesh.vertices:
        if(3*abs(P[0])/np.sqrt(3)>=P[1]):
            ExplicitSolution.append(2*P[1])
        elif(abs(P[0]) < 0.0001):
            ExplicitSolution.append(P[1])
        else:
            ExplicitSolution.append(P[1] + 3*abs(P[0])/np.sqrt(3))

    MySol = Mesh2D.VertexValueMap(MyMesh)
    MySol.values = ExplicitSolution

    Opt=dict(method = "FMT", constrained = False, NarrowBandDepth = 2)

    FMTuSolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt)
    FMTuSolv.computeField()

    Opt3=dict(method = "FMT", constrained = True, NarrowBandDepth = 2)

    FMTCSolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt3)
    FMTCSolv.computeField()

    Opt2=dict(method="FME")

    FMESolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt2)
    FMESolv.computeField()

    for i,VAL in enumerate([FMESolv.fieldValues.values,FMTuSolv.fieldValues.values,FMTCSolv.fieldValues.values]):
        Diff = Mesh2D.VertexValueMap(MyMesh)
        Diff.values = np.array(VAL) - np.array(MySol.values)
        MaxDiffs[i].append(max(Diff.values) if max(Diff.values) > -min(Diff.values) else min(Diff.values))
        Means[i].append(sum(abs(Diff.values))/len(Diff.values))

    for i,VAL in enumerate([FMESolv.fieldValues,FMTuSolv.fieldValues,FMTCSolv.fieldValues]):
        DiffAngles = []
        for j,v in enumerate(VAL.computeGradientFlow()):
            DiffAngles.append(np.acos(v[0]*ConstantVectorField[j][0] + v[1]*ConstantVectorField[j][1]))
        VectorDiffMax[i].append(max(DiffAngles))
        DiffAngles = np.array(DiffAngles)
        VectorDiffMeans[i].append(sum(abs(DiffAngles))/len(DiffAngles))

        DiffAngles = []
        for j,v in enumerate(VAL.computeVertexGradientFlow()):
            DiffAngles.append(abs(np.acos(v[0]*ConstantVectorField[j][0] + v[1]*ConstantVectorField[j][1])))
        VectorDiffMaxE[i].append(max(DiffAngles))
        DiffAngles = np.array(DiffAngles)
        VectorDiffMeansE[i].append(sum(abs(DiffAngles))/len(DiffAngles))

ListPlots = [MaxDiffs,Means,VectorDiffMax,VectorDiffMeans,VectorDiffMaxE,VectorDiffMeansE]

ListNames = ["Maximal signed difference", "Mean difference", "Max angular difference of the gradient", "Mean angular difference of the gradient","Max angular difference of the gradient", "Mean angular difference of the gradient"]

for numPlot in range(6):
    fig, axs = plt.subplots(1)
    for numMod, model in enumerate(['FME', 'FMTU', 'FMTC']):
            axs.plot(X,ListPlots[numPlot][numMod], label = model)
            axs.set_title(ListNames[numPlot] + " with the explicit solution.")
            axs.set_xlabel("1/|Δ|")
            axs.legend()
    #plt.show()
    plt.savefig("figs/plotID2_"+str(numPlot)+".png")
