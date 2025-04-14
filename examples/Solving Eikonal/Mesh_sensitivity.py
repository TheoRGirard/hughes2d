from hughes2d import *
import random as alea
import numpy as np
import matplotlib.pyplot as plt

MyDomain = Mesh2D.NonConvexDomain()
MyDomain = Mesh2D.NonConvexDomain([[0,0],[10,0],[10,5],[0,5]])
MyDomain.addExit([[10,0],[10,5]])

MaxDiffs = [[],[],[],[],[],[]] #order TVE vs explicit, TVTuD2 vs explicit, TVTC vs explicit
Means = [[],[],[],[],[],[]]


X = np.arange(5,50, 5)

for d in X:
    print("Step : ", 1/d)
    MyMesh = Mesh2D.Mesh()
    MyMesh.generateMeshFromDomain(MyDomain, 1/d)
    #MyMesh.show()

    #Custom mesh generate:
    N = 2*int(5*np.sqrt(d)/np.sqrt(2))
    Delta = 10/N
    print("Step : ", Delta**2/2)
    print("Delta : ", Delta)
    verts = []
    tris = []
    Nx = N + 1
    Ny = N//2 + 1

    for n, xn in enumerate([Delta*i for i in range(Nx)]):
        if n%2 == 0:
            for m, yn in enumerate([Delta*j for j in range(Ny)]):
                verts.append([xn,yn])
                if(n>0):
                    FirstLastLine = (n-1)*Ny+(n-1)//2
                    FirstThisLine = n*Ny+n//2
                    tris.append([FirstLastLine +m,FirstLastLine +m+1,FirstThisLine +m])
                    if(m < Ny-1):
                        tris.append([FirstLastLine +m+1,FirstThisLine +m,FirstThisLine+m+1])
        else:
            for m,yn in enumerate([0]+[Delta/2 + Delta*j for j in range(Ny-1)] + [5]):
                verts.append([xn,yn])
                if(m < Ny):
                    FirstLastLine = (n-1)*Ny+(n-1)//2
                    FirstThisLine = n*Ny+n//2
                    if(m > 0 and m < N+1):
                        tris.append([FirstLastLine+m-1,FirstLastLine+m,FirstThisLine+m])
                    tris.append([FirstLastLine+m,FirstThisLine+m,FirstThisLine+m+1])

    MyMesh2 = Mesh2D.Mesh()
    print("Nb triangles : ", len(tris))

    MyMesh2.importFromLists(vertices = verts,triangles = tris,domain = MyDomain)

    #MyMesh2.show()

    MyMap = Mesh2D.CellValueMap(MyMesh)
    MyMap2 = Mesh2D.CellValueMap(MyMesh2)

    ConstantVectorField = []
    ConstantVectorField2 = []

    for P in MyMesh.triangles:
        ConstantVectorField.append([1,0])
    for P in MyMesh2.triangles:
        ConstantVectorField2.append([1,0])

    ExplicitSolution = []
    ExplicitSolution2 = []

    for P in MyMesh.vertices:
        ExplicitSolution.append(10-P[0])
    for P in MyMesh2.vertices:
        ExplicitSolution2.append(10-P[0])

    MySol = Mesh2D.VertexValueMap(MyMesh)
    MySol.values = ExplicitSolution
    MySol2 = Mesh2D.VertexValueMap(MyMesh2)
    MySol2.values = ExplicitSolution2

    Opt=dict(method = "FMT", constrained = False, NarrowBandDepth = 2)

    FMTuSolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt)
    FMTuSolv.computeField()
    FMTuSolv2 = EikonalSolver.EikoSolver(MyMesh2, DensityMap = MyMap2, opt = Opt)
    FMTuSolv2.computeField()


    Opt3=dict(method = "FMT", constrained = True, NarrowBandDepth = 2)

    FMTCSolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt3)
    FMTCSolv.computeField()
    FMTCSolv2 = EikonalSolver.EikoSolver(MyMesh2, DensityMap = MyMap2, opt = Opt3)
    FMTCSolv2.computeField()

    Opt2=dict(method="FME")

    FMESolv = EikonalSolver.EikoSolver(MyMesh, DensityMap = MyMap, opt = Opt2)
    FMESolv.computeField()
    FMESolv2 = EikonalSolver.EikoSolver(MyMesh2, DensityMap = MyMap2, opt = Opt2)
    FMESolv2.computeField()

    for i,VAL in enumerate([FMESolv.fieldValues.values,FMTuSolv.fieldValues.values,FMTCSolv.fieldValues.values]):
        Diff = Mesh2D.VertexValueMap(MyMesh)
        Diff.values = np.array(VAL) - np.array(MySol.values)
        MaxDiffs[i].append(max(Diff.values) if max(Diff.values) > -min(Diff.values) else min(Diff.values))
        Means[i].append(sum(abs(Diff.values))/len(Diff.values))

    for i,VAL in enumerate([FMESolv2.fieldValues.values,FMTuSolv2.fieldValues.values,FMTCSolv2.fieldValues.values]):
        Diff = Mesh2D.VertexValueMap(MyMesh2)
        Diff.values = np.array(VAL) - np.array(MySol2.values)
        MaxDiffs[i+3].append(max(Diff.values) if max(Diff.values) > -min(Diff.values) else min(Diff.values))
        Means[i+3].append(sum(abs(Diff.values))/len(Diff.values))

ListPlots = [MaxDiffs,Means]

ListNames = ["Maximal signed difference", "Mean difference"]
colors = ['b','C1','g']

for numPlot in range(2):
    fig, axs = plt.subplots(1)
    for numMod, model in enumerate(['FME-regular mesh', 'FMTU-regular mesh', 'FMTC-regular mesh']):
            axs.plot(X,ListPlots[numPlot][numMod+3],"-", label = model)

    for numMod, model in enumerate(['FME', 'FMTU', 'FMTC']):
            axs.plot(X,ListPlots[numPlot][numMod],colors[numMod]+ "--", label = model)

    axs.set_title(ListNames[numPlot] + " with the explicit solution.")
    axs.set_xlabel("1/|Δ|")
    axs.legend()
    #plt.show()
    plt.savefig("figs/plotInfluenceMesh"+str(numPlot)+".png")
