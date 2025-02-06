import numpy as np
import numpy.typing as npt

Point = np.array([0.0,1.2])

Triangle = np.array([1,0,3])

ListOfPoint = np.array([[0.1,5],[4.1,8],[6,0.2]])

print("Type of Point : ", type(Point))
print("Type of ListOfPoint : ", type(ListOfPoint))

print(npt.NDArray)
print(npt.NDArray[int])
print(npt.NDArray[float])
