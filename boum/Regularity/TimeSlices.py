from hughes2d.Plotter import *

saveTimeSlices([0,1,2]+ [3+i*3 for i in range(9)],filename="data/TestRegu_1",slicename = "figs/HughesRegu_", limits=[[0,10],[0,5]])
