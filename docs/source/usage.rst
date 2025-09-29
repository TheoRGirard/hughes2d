Usage
=====

.. _installation:

Installation
------------

The hughes2d package dependencies are managed using uv (see uv). If uv is installed, you can install the hughes2d package and its dependencies by typing:

.. code-block:: console

   git clone https://github.com/TheoRGirard/hughes2d |
   cd hughes2d |
   uv sync

If you don't have uv installed, you can see a list of the dependencies in the pyproject.toml file or below.

hughes2d v0.1.0
├── numpy v2.2.2
├── triangle v20250106
│   └── numpy v2.2.2
├── ezdxf v1.3.5 (extra: file-formats)
│   ├── fonttools v4.55.8
│   ├── numpy v2.2.2
│   ├── pyparsing v3.2.1
│   └── typing-extensions v4.12.2
├── meshio v5.3.5 (extra: file-formats)
│   ├── numpy v2.2.2
│   └── rich v13.9.4
│       ├── markdown-it-py v3.0.0
│       │   └── mdurl v0.1.2
│       └── pygments v2.19.1
├── matplotlib v3.10.0 (extra: plot)
│   ├── contourpy v1.3.1
│   │   └── numpy v2.2.2
│   ├── cycler v0.12.1
│   ├── fonttools v4.55.8
│   ├── kiwisolver v1.4.8
│   ├── numpy v2.2.2
│   ├── packaging v24.2
│   ├── pillow v11.1.0
│   ├── pyparsing v3.2.1
│   └── python-dateutil v2.9.0.post0
│       └── six v1.17.0
├── plotly v6.0.0 (extra: plot)
│   ├── narwhals v1.24.1
│   └── packaging v24.2
├── pyqt6 v6.8.1 (extra: plot)
│   ├── pyqt6-qt6 v6.8.2
│   └── pyqt6-sip v13.10.0



Getting Started
----------------

You can find a file named getting_started.py in /examples. We rewrite below the content of this file.

::

   from hughes2d import *
   
   #Construction of the domain--------------------------------
   MyDomain = Mesh2D.NonConvexDomain([[0,0],[0,1],[1,1],[1,0]])
   MyDomain.addExits([[[1,0],[1,1]]])
   MyDomain.show()
   
   #Construction of the mesh--------------------------------------
   MyMesh = Mesh2D.Mesh()
   MyMesh.generateMeshFromDomain(MyDomain, 0.01)
   MyMesh.show()
   MyMesh.saveToJson("gettingStartedSimu")
   
   
   #Construction of a random initial datum---------------------------------------
   MyMap = Mesh2D.CellValueMap(MyMesh)
   MyMap.generateRandom()
   MyMap.show()
   
   #Setting the options for the simulation-----------------------------------------
   opt = dict(model = "hughes",
               filename = "gettingStartedSimu",
               save = True,
               verbose = True
               )
   
   #Creating the solver and computing---------------------------------------------------
   Solver = Splitting.PedestrianSolver(MyMesh, 0.01,0.01, initialDensity = MyMap, options=opt)
   Solver.computeUntilEmpty(100)
   
   #Converting the data to a mp4 video------------------------------------------
   Plotter.convertToMP4("gettingStartedSimu", limits=[[0,1],[0,1]])


Compiling and running this code should create 3 .csv files, 1 .json file and 1 .mp4 file. The .mp4 file should look like the video below:

.. image: gettingStartedVid.gif
