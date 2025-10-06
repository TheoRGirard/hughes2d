
"""
Solver for LWR-like 2D equations

Girard Theo
"""


from hughes2d.Mesh2D import *

import copy
import json

class LWRSolver(object):
    """
    A solver for LWR-type scalar conservation law directed by a vector field in two dimensions. We use the following notations:

    .. math::

        \\left\\{ \\begin{matrix}
        \\partial_t \\rho(t,x) + \\mathbf{div}\\left[ \\rho(t,x) v(\\rho(t,x)) \\vec{V}(t,x) \\right] = 0, \\\\
        \\rho(0,x) = \\rho_0(x)
        \\end{matrix} \\right.

    Args:
        Mesh (Mesh): a mesh object on which the equation will be approximated.
        previousDensity (CellValueMap or List[float]): initial density for the solver. Must be of the shape of `Mesh.triangles`. Represents :math:`\\rho_0(x)` in the equation above.
        directionMap (List[List[float]]): a vector field represented by a list of vectors with the shape of `Mesh.triangles`. Typically this corresponds to the output of `VertexValueMap.computeGradientFlow()`. Represents :math:`\\vec{V}(t,x)` in the equation above.
        speedFunction (function, float -> float): a function respresenting the speed of the agents depending on the local density. Represents :math:`v(\\cdot)` in the equation above.
        dt (float): the time division for the approximation.
            Warning:
                The CFL condition must be satisfied for the simulations to make sense. Here the CFL condition is:

                .. math::
                    \\Delta t \\leq \\frac{\\underline{|\\triangle|}}{3\\underline{\\textrm{$\\triangle$}}Lip_f},

                where :math:`\\underline{|\\triangle|}` denotes the minimal area of a triangle in the mesh :math:`M_\\Delta` and :math:`\\underline{\\textrm{$\\triangle$}}` denotes the maximal length of the edges of the mesh.
        opt (dict): an optional dictionary prescribing the numerical method.

            ========================= ====== ====================== =========================================================================================
            key                       type   possible values        description
            ========================= ====== ====================== =========================================================================================
            'method'                  str    "tmap" or "midvector"  determines of the conflict between non-colinear vectors are resolved
            'anNum'                   str    "dichotomy"            only parameter available at the moment, numerical method to use for the approximations
            'convexFlux'              bool   True or False          optimization of the computations when the flux is convex or concave
            'ApproximationThreshold'  float  > 0                    the precision to use for the numerical approximations
            'debugging'               bool   True or False          determines if the solver will print debugging informations in the console or not
            ========================= ====== ====================== =========================================================================================

    Attributes
    -----------

    Attributes:
        mesh (Mesh): a mesh object on which the equation will be approximated.
        densityt0 (ArrayLike): initial density for the solver. Of the shape of `Mesh.triangles`. Represents :math:`\\rho_0(x)` in the equation above.
        densityt1 (ArrayLike): density computed by the solver after one time step. Of the shape of `Mesh.triangles`. Represents :math:`\\rho(\\Delta t, x)`.
        directions (List[List[float]]): a vector field represented by a list of vectors with the shape of `Mesh.triangles`. Typically this corresponds to the output of `VertexValueMap.computeGradientFlow()`. Represents :math:`\\vec{V}(t,x)` in the equation above.
        fluxFunction (function, float -> float): a function respresenting the flux of agents depending on the local density. Represents :math:`\\rho \\mapsto \\rho v(\\cdot)` in the equation above.
        dt (float): the time division for the approximation.

            Warning:
                The CFL condition must be satisfied for the simulations to make sense. Here the CFL condition is:

                .. math::
                    \\Delta t \\leq \\frac{\\underline{|\\triangle|}}{3\\underline{\\textrm{$\\triangle$}}Lip_f},

                where :math:`\\underline{|\\triangle|}` denotes the minimal area of a triangle in the mesh :math:`M_\\Delta` and :math:`\\underline{\\textrm{$\\triangle$}}` denotes the maximal length of the edges of the mesh.

        opt (dict): an optional dictionary prescribing the numerical method.

    Methods
    -----------
    """

    def __init__(self, Mesh:Mesh, dt:float, previousDensity=[], directionMap:List[List[float]] = [], speedFunction = (lambda x: 1-x), opt:dict = dict(convexFlux = True, anNum = "dichotomy", method = "midVector")):
        self.mesh:Mesh = Mesh
        if(previousDensity != []):
            self.densityt0:ArrayLike = np.array(previousDensity) #type CellValueMap
        else:
            self.densityt0:ArrayLike = np.array([0.0 for t in self.mesh.triangles])
        self.densityt1:ArrayLike = np.array([0.0 for t in self.mesh.triangles])

        self.dt:float = dt

        self.opt:dict = opt

        self.directions:List[List[float]] = directionMap

        if("ApproximationThreshold" not in self.opt.keys() ):
            self.opt['ApproximationThreshold'] = 0.0001

        if("debugging" not in self.opt.keys() ):
            self.opt['debugging'] = False

        self.fluxFunction = lambda x: x*speedFunction(x)

        if(self.opt["convexFlux"]):
            if(self.fluxFunction == (lambda x: x*(1-x))):
                self.maxFluxPoint = 0.5
                if self.opt['debugging']:
                    print("Explicit max point considered to be 0.5")
            else:
                self.maxFluxPoint = argMax(self.fluxFunction,0,1,self.opt['ApproximationThreshold'])
                if self.opt['debugging']:
                    print("Maximal flux point found at p = ", self.maxFluxPoint)

    def checkCFL(self, LipConstant:float = 1) -> bool:
        """
        Checks if the CFL condition is satisfied. Here the CFL condition is:

        .. math::
            \\Delta t \\leq \\frac{\\underline{|\\triangle|}}{3\\underline{\\textrm{$\\triangle$}}Lip_f},

        where :math:`\\underline{|\\triangle|}` denotes the minimal area of a triangle in the mesh :math:`M_\\Delta` and :math:`\\underline{\\textrm{$\\triangle$}}` denotes the maximal length of the edges of the mesh.

        Args:
            LipConstant (float): the Lipschitz constant of the flux.

        Returns:
            bool: True if the CFL condition is verified.
        """
        if 3*self.dt*LipConstant*self.mesh.maxEdgeLength > self.mesh.minCellArea:
            print("Warning - the CFL conditions might not be verified everywhere in the mesh.")
            if self.opt['debugging']:
                print("CFL values : dt = ", self.dt, " VS ", self.mesh.minCellArea/(3*LipConstant*self.mesh.maxEdgeLength))
            return False
        return True


    def computeStepTmap(self) -> None:
        """
        Computes the approximated solution to the scalar conservation law after one time step with the transmission maps method.
        """
        if(self.opt["convexFlux"]):
            for triangleIndex, triangleCell in enumerate(self.mesh.trianglesWithEdges):
                Modifdensity = 0
                for edgeNumber, edgeIndex in enumerate(triangleCell):
                    PairOfTriangles = self.mesh.pairsOfTriangles[edgeIndex]
                    TriangleGrad = self.directions[triangleIndex]
                    Normal = self.mesh.outerNormalVectByTriangles[triangleIndex][edgeNumber]

                    triangleFlux = TriangleGrad[0]*Normal[0] + TriangleGrad[1]*Normal[1]

                    if(edgeIndex in self.mesh.exitEdges):
                        farFlux = 1
                        farDensity = 0
                    elif(edgeIndex in self.mesh.wallEdges):
                        triangleFlux = 0
                        farFlux = 0
                        farDensity = 0
                    else:
                        if(PairOfTriangles[0] == triangleIndex):
                            farTriangleGrad = self.directions[PairOfTriangles[1]]
                            farDensity = self.densityt0[PairOfTriangles[1]]
                        else:
                            farTriangleGrad = self.directions[PairOfTriangles[0]]
                            farDensity = self.densityt0[PairOfTriangles[0]]
                        farFlux = (farTriangleGrad[0]*Normal[0] + farTriangleGrad[1]*Normal[1])

                    if (triangleFlux > 1 + self.opt['ApproximationThreshold'] or farFlux > 1 + self.opt['ApproximationThreshold']):
                        print("Warning : scalar product out of bonds...")
                        if self.opt['debugging']:
                            print("| farTriangleGrad | = ", farTriangleGrad[0]*farTriangleGrad[0] + farTriangleGrad[1]*farTriangleGrad[1])
                            print("| normal | = ", Normal[0]*Normal[0] + Normal[1]*Normal[1])


                    """
                    We search for k such that God(triangleFlux function, densityt0[triangleIndex], k ) = God(farFlux function, k, fardensity).
                    In the specific case of a convex flux function we can solve that by treating different cases :
                    """
                    if(triangleFlux <= 0 and farFlux >= 0):
                        Totalflux = 0
                    elif(triangleFlux >= 0 and farFlux <= 0):
                        Totalflux = 0

                    elif(triangleFlux < 0 and farFlux < 0):
                        OutingFlux = triangleFlux*self.fluxFunction(max(self.densityt0[triangleIndex],self.maxFluxPoint))
                        EnteringFlux = farFlux*self.fluxFunction(min(farDensity,self.maxFluxPoint))
                        Totalflux = max(OutingFlux,EnteringFlux)

                    elif(triangleFlux > 0 and farFlux > 0):
                        EnteringFlux = farFlux*self.fluxFunction(max(farDensity,self.maxFluxPoint))
                        OutingFlux = triangleFlux*self.fluxFunction(min(self.densityt0[triangleIndex],self.maxFluxPoint))
                        Totalflux = min(OutingFlux,EnteringFlux)
                    Modifdensity -= (self.dt/self.mesh.cellAreas[triangleIndex]) * self.mesh.edgeLength[edgeIndex]* Totalflux
                self.densityt1[triangleIndex] = self.densityt0[triangleIndex] + Modifdensity

        else: #flux not convex
            for triangleIndex, triangleCell in enumerate(self.mesh.trianglesWithEdges):
                Modifdensity = 0
                for edgeNumber, edgeIndex in enumerate(triangleCell):
                    PairOfTriangles = self.mesh.pairsOfTriangles[edgeIndex]
                    TriangleGrad = self.directions[triangleIndex]
                    Normal = self.mesh.outerNormalVectByTriangles[triangleIndex][edgeNumber]

                    triangleFlux = TriangleGrad[0]*Normal[0] + TriangleGrad[1]*Normal[1]

                    if(edgeIndex in self.mesh.exitEdges):
                        farFlux = 1
                        farDensity = 0
                    elif(edgeIndex in self.mesh.wallEdges):
                        triangleFlux = 0
                        farFlux = 0
                        farDensity = 0
                    else:
                        if(PairOfTriangles[0] == triangleIndex):
                            farTriangleGrad = self.directions[PairOfTriangles[1]]
                            farDensity = self.densityt0[PairOfTriangles[1]]
                        else:
                            farTriangleGrad = self.directions[PairOfTriangles[0]]
                            farDensity = self.densityt0[PairOfTriangles[0]]
                        farFlux = (farTriangleGrad[0]*Normal[0] + farTriangleGrad[1]*Normal[1])

                    if(triangleFlux > 1 + self.opt['ApproximationThreshold'] or farFlux > 1 + self.opt['ApproximationThreshold']):
                        print("Warning : scalar product out of bounds")
                        if self.opt['debugging']:
                            print("| farTriangleGrad | = ", farTriangleGrad[0]*farTriangleGrad[0] + farTriangleGrad[1]*farTriangleGrad[1])
                            print("| normal | = ", Normal[0]*Normal[0] + Normal[1]*Normal[1])
                    """
                    If the flux function is not convex, we use a dichotomy to solve the problem.
                    """
                    if(triangleFlux <= 0 and farFlux >= 0):
                        Totalflux = 0
                    elif(triangleFlux >= 0 and farFlux <= 0):
                        Totalflux = 0
                    else:
                        Outerflux = lambda x: triangleFlux*self.fluxFunction(x)
                        Innerflux = lambda x: farFlux*self.fluxFunction(x)
                        parametrizedFlux = lambda k : God(Outerflux, self.densityt0[triangleIndex], k,precision=self.opt['ApproximationThreshold']) - God(Innerflux,k, farDensity, precision=self.opt['ApproximationThreshold'])
                        k = ApproZeroDichotomie(parametrizedFlux,0,1,self.opt['ApproximationThreshold'], hints=[self.densityt0[triangleIndex],farDensity])
                        Totalflux = God(Outerflux, self.densityt0[triangleIndex], k,self.opt['ApproximationThreshold'])
                    Modifdensity -= (self.dt/self.mesh.cellAreas[triangleIndex]) * self.mesh.edgeLength[edgeIndex]* Totalflux

                self.densityt1[triangleIndex] = min(1,max(self.densityt0[triangleIndex] + Modifdensity,0))

                if(np.abs(min(1,max(self.densityt0[triangleIndex] + Modifdensity,0)) - (self.densityt0[triangleIndex] + Modifdensity)) > float(1e-10)):
                    print("Warning : the new density computed is not between 0 and 1.")
                    if(self.opt['debugging']):
                        print("Computed density: ", self.densityt0[triangleIndex] + Modifdensity)

    def computeStepMidVector(self) -> None:
        """
        Computes the approximated solution to the scalar conservation law after one time step with the mid vector method.
        """
        if(self.opt["convexFlux"]):
            for triangleIndex, triangleCell in enumerate(self.mesh.trianglesWithEdges):
                Modifdensity = 0

                for edgeNumber, edgeIndex in enumerate(triangleCell):

                    PairOfTriangles = self.mesh.pairsOfTriangles[edgeIndex]
                    TriangleGrad = self.directions[triangleIndex]
                    Normal = self.mesh.outerNormalVectByTriangles[triangleIndex][edgeNumber]
                    Totalflux = 1

                    otherTriangleIndex = -1

                    if(edgeIndex in self.mesh.exitEdges):
                        VectorFlux = TriangleGrad
                        farDensity = 0
                    elif(edgeIndex in self.mesh.wallEdges):
                        VectorFlux = TriangleGrad
                        Totalflux = 0
                    else:

                        if(PairOfTriangles[0] == triangleIndex):
                            farTriangleGrad = self.directions[PairOfTriangles[1]]
                            farDensity = self.densityt0[PairOfTriangles[1]]

                            otherTriangleIndex = PairOfTriangles[1]
                        else:
                            farTriangleGrad = self.directions[PairOfTriangles[0]]
                            farDensity = self.densityt0[PairOfTriangles[0]]

                            otherTriangleIndex = PairOfTriangles[0]

                        VectorFlux = [(TriangleGrad[0]*self.densityt0[triangleIndex] + farTriangleGrad[0]*farDensity), (TriangleGrad[1]*self.densityt0[triangleIndex] + farTriangleGrad[1]*farDensity) ]
                    normeVectorFlux = np.sqrt(VectorFlux[0]*VectorFlux[0] + VectorFlux[1]*VectorFlux[1])

                    if(normeVectorFlux < self.opt['ApproximationThreshold'] or Totalflux == 0):
                        Totalflux = 0
                    else:
                        Scalar = (VectorFlux[0]*Normal[0] + VectorFlux[1]*Normal[1])
                        if(Scalar/normeVectorFlux > 1 + float(1e-10)):
                            print("Warning : scalar product greater than the norm, perhaps the normal vector is not normalized")
                        if(Scalar > 0):
                            if(self.densityt0[triangleIndex] <= farDensity):
                                Totalflux = Scalar*min(self.fluxFunction(self.densityt0[triangleIndex]),self.fluxFunction(farDensity))/normeVectorFlux
                            elif(farDensity < self.maxFluxPoint and self.maxFluxPoint < self.densityt0[triangleIndex]):
                                Totalflux = Scalar*self.fluxFunction(self.maxFluxPoint)/normeVectorFlux
                            else:
                                Totalflux = Scalar*max(self.fluxFunction(self.densityt0[triangleIndex]),self.fluxFunction(farDensity))/normeVectorFlux
                        else:
                            if(self.densityt0[triangleIndex] >= farDensity):
                                Totalflux = Scalar*min(self.fluxFunction(self.densityt0[triangleIndex]),self.fluxFunction(farDensity))/normeVectorFlux
                            elif(farDensity > self.maxFluxPoint and self.maxFluxPoint > self.densityt0[triangleIndex]):
                                Totalflux = Scalar*self.fluxFunction(self.maxFluxPoint)/normeVectorFlux
                            else:
                                Totalflux = Scalar*max(self.fluxFunction(self.densityt0[triangleIndex]),self.fluxFunction(farDensity))/normeVectorFlux

                    Modifdensity -= (self.dt/self.mesh.cellAreas[triangleIndex]) * self.mesh.edgeLength[edgeIndex]* Totalflux

                self.densityt1[triangleIndex] = min(1,max(self.densityt0[triangleIndex] + Modifdensity,0))

                if(np.abs(min(1,max(self.densityt0[triangleIndex] + Modifdensity,0)) - (self.densityt0[triangleIndex] + Modifdensity)) > float(1e-10)):
                    print("Warning : the new density computed is not between 0 and 1.")
                    if(self.opt['debugging']):
                        print("Computed density: ", self.densityt0[triangleIndex] + Modifdensity)

        else: #for non convex flux
            for triangleIndex, triangleCell in enumerate(self.mesh.trianglesWithEdges):
                Modifdensity = 0
                for edgeNumber, edgeIndex in enumerate(triangleCell):
                    PairOfTriangles = self.mesh.pairsOfTriangles[edgeIndex]
                    TriangleGrad = self.directions[triangleIndex]
                    Normal = self.mesh.outerNormalVectByTriangles[triangleIndex][edgeNumber]
                    Totalflux = 1

                    if(edgeIndex in self.mesh.exitEdges):
                        VectorFlux = TriangleGrad
                        farDensity = 0
                    elif(edgeIndex in self.mesh.wallEdges):
                        VectorFlux = TriangleGrad
                        Totalflux = 0
                    else:

                        if(PairOfTriangles[0] == triangleIndex):
                            farTriangleGrad = self.directions[PairOfTriangles[1]]
                            farDensity = self.densityt0[PairOfTriangles[1]]
                        else:
                            farTriangleGrad = self.directions[PairOfTriangles[0]]
                            farDensity = self.densityt0[PairOfTriangles[0]]

                        VectorFlux = [(TriangleGrad[0]*self.densityt0[triangleIndex] + farTriangleGrad[0]*farDensity), (TriangleGrad[1]*self.densityt0[triangleIndex] + farTriangleGrad[1]*farDensity) ]
                    normeVectorFlux = np.sqrt(VectorFlux[0]*VectorFlux[0] + VectorFlux[1]*VectorFlux[1])

                    if(normeVectorFlux < self.opt['ApproximationThreshold'] or Totalflux == 0):
                        Totalflux = 0
                    else:
                        if((VectorFlux[0]*Normal[0] + VectorFlux[1]*Normal[1])/normeVectorFlux > 1 + float(1e-10)):
                            print("Warning : scalar product greater than the norm, perhaps the normal vector is not normalized")
                        ModifFluxFunc = lambda x : (VectorFlux[0]*Normal[0] + VectorFlux[1]*Normal[1])/normeVectorFlux *self.fluxFunction(x)
                        Totalflux = God(ModifFluxFunc, self.densityt0[triangleIndex], farDensity, precision=self.opt['ApproximationThreshold'])

                    Modifdensity -= (self.dt/self.mesh.cellAreas[triangleIndex]) * self.mesh.edgeLength[edgeIndex]* Totalflux
                self.densityt1[triangleIndex] = min(1,max(self.densityt0[triangleIndex] + Modifdensity,0))

                if(np.abs(min(1,max(self.densityt0[triangleIndex] + Modifdensity,0)) - (self.densityt0[triangleIndex] + Modifdensity)) > float(1e-10)):
                    print("Warning : the new density computed is not between 0 and 1.")
                    if(self.opt['debugging']):
                        print("Computed density: ", self.densityt0[triangleIndex] + Modifdensity)

    def computeNextStep(self):
        """
        Computes the approximated solution to the scalar conservation law after one time step and stores it in `densityt1`.
        """
        if(self.opt["method"] == "tmap"):
            return self.computeStepTmap()
        elif(self.opt["method"] == "midVector"):
            return self.computeStepMidVector()

    def update(self, newDirectionField:List[List[float]]) -> None:
        """
        Updates the direction vector field with the one passed as a parameter and sets the actual approximated solution `densityt1` as the initial datum `densityt0`.

        Args:
            newDirectionField (List[List[float]]): a vector field represented by a list of vectors with the shape of `Mesh.triangles`. Typically this corresponds to the output of `VertexValueMap.computeGradientFlow()`. Represents :math:`\\vec{V}(t,x)` in the scalar conservation law equation.
        """
        self.densityt0 = np.copy(self.densityt1)
        self.directions = newDirectionField

def argMax(f,a:float,b:float, precision:float = 0.0001):
    """
    Numerically approximates the maximal point of the function `f` between `a` and `b` and returns the maximal argument.

    Args:
        f (function, float -> float): the function `f` for which the armax will be computed.
        a (float): one of the bounds for the search domain of the argmax.
        b (float): one of the bounds for the search domain of the argmax.
        precision (float): the error margin for the approximation.

    Returns:
        float: the argmax computed.
    """
    NumSlice = int(1+1/precision)
    Max = -float('inf')
    Pas = abs(b-a)/NumSlice
    xmax = min(a,b)

    for i in range(NumSlice):
        Test = f(min(a,b) + i*Pas)
        if(Test > Max):
            Max = Test
            xmax = min(a,b) + i*Pas
    return(xmax)

def Max(f,a:float,b:float, precision:float = 0.0001):
    """
    Numerically approximates the maximum of the function `f` between `a` and `b` and returns the maximal value.

    Args:
        f (function, float -> float): the function `f` for which the max will be computed.
        a (float): one of the bounds for the search domain of the max.
        b (float): one of the bounds for the search domain of the max.
        precision (float): the error margin for the approximation.

    Returns:
        float: the maximal value computed.
    """
    Max = -float('inf')
    nbStep:int = int(np.ceil(np.abs(b-a)/precision))
    if nbStep == 0:
        print("Warning: precision is less than the difference between the bounds of the domain")
        return (a+b)/2

    Pas = np.abs(b-a)/nbStep

    for i in range(nbStep):
        Test = f(min(a,b) + i*Pas)
        if(Test > Max):
            Max = Test
    return(Max)

def Min(f,a:float,b:float, precision:float = 0.0001):
    """
    Numerically approximates the minimum of the function `f` between `a` and `b` and returns the minimal value.

    Args:
        f (function, float -> float): the function `f` for which the min will be computed.
        a (float): one of the bounds for the search domain of the min.
        b (float): one of the bounds for the search domain of the min.
        precision (float): the error margin for the approximation.

    Returns:
        float: the minimal value computed.
    """
    Min = float('inf')
    nbStep:int = int(np.ceil(np.abs(b-a)/precision))
    if nbStep == 0:
        print("Warning: precision is less than the difference between the bounds of the domain")
        return (a+b)/2

    Pas = np.abs(b-a)/nbStep

    for i in range(nbStep):
        Test = f(min(a,b) + i*Pas)
        if(Test < Min):
            Min = Test
    return(Min)

def God(f,a:float,b:float, precision:float = 0.0001):
    """
    Numerically approximates the Godunov flux of the function `f` between `a` and `b` defined by the formula:

    .. math::

        \\mathbf{God}_{f}(a,b) = \\left\\{ \\begin{matrix} \\mathbf{min}_{c \\in [a,b]} f(c) \\textrm{ if } a < b \\\\
        \\mathbf{max}_{c \\in [b,a]} f(c) \\textrm{ if } a > b. \\end{matrix} \\right.

    Args:
        f (function, float -> float): the function `f` for which the Godunov flux will be computed.
        a (float): one of the parameters for the Godunov flux.
        b (float): one of the parameters for the Godunov flux.
        precision (float): the error margin for the approximation.

    Returns:
        float: the Godunov flux computed.
    """
    if np.abs(a - b) < precision:
        return f((a+b)/2)
    if(a < b):
        return Min(f,a,b, precision)
    else :
        return Max(f,b,a, precision)

def ApproZeroDichotomie(f,a:float,b:float, precision:float = 0.0001, hints=[]):
    """
    Numerically approximates the a root of the function `f` between `a` and `b` using a dichotomy method.

    Args:
        f (function, float -> float): the function `f` for which a root will be computed.
        a (float): one of the bounds for the search domain of the root.
        b (float): one of the bounds for the search domain of the root.
        precision (float,optional): the error margin for the approximation.
        hints (float,optional): possibles value to test before the dichotomy in order to optimize the computation time.

    Returns:
        float: the approximated root computed.
    """
    for x in hints:
        if abs(f(x)) < precision:
            return x

    if abs(f(a)) < precision:
        return a
    if abs(f(b)) < precision:
        return b

    c = a+ (b-a)/2
    while b-a > precision:
        if abs(f(c)) < precision:
            return c
        if(f(a) > 0 and f(b) < -0):
            if(f(c) >= 0):
                a = c
            else:
                b = c
        elif(f(b) > 0 and f(a) < -0):
            if(f(c)>0):
                b = c
            else:
                a = c
        else:
            raise ValueError("Dichotomy is impossible with the function and the domain given as parameters. You should check if the function change of sign exactly once in the domain.")

        c = a+(b-a)/2

    if abs(f(c)) < abs(f(a)) and abs(f(c)) < abs(f(b)):
        return c
    elif abs(f(a)) < abs(f(b)):
        return a
    else:
        return b
