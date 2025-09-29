Mathematics
=====

.. _macroPedModels:

Macroscopic pedestrian models
-----------------------------

The mathematical models for pedestrian crowds can be catyegorized in three main categories.
- The **microscopic** models use a huge number of agents interacting with each other to model a crowd of pedestrian. This type of model often relies on the use of a system coupling many Ordinary Differential Equation (ODE).

- The **mesoscopic** models deal with an infinite number of agents characterized by an artificial parameter. This type of model mainly consists in kinetic partial differential equations and uses the formalism of statistic physics.

- The **macroscopic** models deal with a density of pedestrian instead of agents. This type of model relies on the use of Partial Differential Equation (PDE) which are often close to the PDEs found in fluid dynamics.

The present python package focuses on the numerical approximation of solutions to Hughes' model (introduced in [Hug02]_ ) which is one of most studied the **macroscopic** models.

Hughes' model
-------------
.. note::
  The introduction to Hughes' model below is quoted from the introduction of Chapter 2 of [Gir25]_. It is recommanded reading it there as more detail are included.

In [Hug02]_, Hughes proposed a mathematical model for the two-dimensional dynamics of pedestrian crowds. The model consists in a system of two equations set on a bounded domain :math:`\Omega \subset \mathbb{R}^2`. The first PDE is a vector-directed *scalar conservation law* with discontinuous flux. The second is an *Eikonal equation* with discontinuous source term.

The scalar conservation law 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first equation models the flow of a density :math:`\rho`. To be more precise  :math:`\rho(t,x)` represents the density of pedestrian at time :math:`t` and location :math:`x`. It is bounded between :math:`0` and a given :math:`\rho_{max} >0`. We assume that the pedestrian move with the speed of agents :math:`v(t,x)` at time :math:`t` and location :math:`x` following a unitary direction field :math:`\vec{V}(t,x) \in \mathcal{S}_1`.
Then, if we write the conservation of the mass on pedestrian on each subdomain of :math:`\Omega`, we end up with the following scalar conservation law:

.. math::
  \partial_t \rho + \mathbf{div}(\vec{V}(t,x) v(t,x) \rho(t,x)) = 0.

In Hughes' model, we assume that the speed of agents :math:`v(t,x)` at time :math:`t` and location :math:`x` only depends on the density :math:`\rho(t,x)` and that this speed is decreasing with respect to the density.
Then, in the following we will denote the speed by :math:`v(\rho(t,x))` where :math:`v` is a decreasing function defined on :math:`[0,\rho_{max}]` such that

.. math:: 
  v(0)=: v_{\max} \textrm{ and } v(\rho_{max})= 0.

A classical example of such a speed function is :math:`v(\rho) = v_{max}\frac{\rho_{max}-\rho}{\rho_{max}}`. This corresponds to the very classical Lighthill-Whitham-Richards (LWR) model for traffic flows (see [LW55]_, [Ric56]_).

The Eikonal equation
^^^^^^^^^^^^^^^^^^^^^

The second equation of the model characterizes the unitary direction field :math:`\vec{V}(t,x)` depending on the density :math:`\rho` in the whole domain. We assume that the pedestrians want to minimize their exit time while also trying to avoid high density regions. In order to model this situation, we use an optimal control problem. We suppose that the density :math:`\rho(\cdot) \in \mathcal{C}^1(\bar \Omega)` stays constant in time (this assumption is quite controversial).
Let :math:`x \in \Omega`. For any :math:`\alpha(\cdot) \in \mathcal{C}^1((0,+\infty),\mathcal{S}_1)`, we say that :math:`X^\alpha_x(\cdot)` is a trajectory controlled by :math:`\alpha` starting at :math:`x` if :math:`X` is a solution to the Cauchy problem:

.. math::
  \left\lbrace \begin{matrix}
  \dot{X}^\alpha_x(s) = \alpha(s)v(\rho(t,X^\alpha_x(s))) \\
  X^\alpha_x(0) = x
  \end{matrix}\right.

For any :math:`x \in \Omega`, we denote by :math:`\mathcal{A}_x = \{ X^\alpha_x, \alpha \in \mathcal{C}^1((0,+\infty),\mathcal{S}_1) \}` the set of all controlled trajectories starting at :math:`x`. We define :math:`\phi(x)` as the minimal exit time starting at location :math:`x`, that is to say:

.. math::
  \phi(x) = \inf_{X \in \mathcal{A}_x} \int_0^{+\infty} \mathbb{1}_{\Omega}(X(t)) \textrm{d} t.

In Hughes' model, we also take into account the discomfort caused by being surrounded by a high density crowd. In order to model this discomfort, we introduce an increasing function :math:`g(\rho)` with respect to the density :math:`\rho`. The function :math:`g(\rho)` can be interpreted as a running cost we are paying along a trajectory :math:`X` for being in high density regions. Then, the previous equation becomes:

.. math::
  \phi(x) = \inf_{X \in \mathcal{A}_x} \int_0^{+\infty} \mathbb{1}_{\Omega}(X(t))g(\rho(X(t))) 	extrm{d} t.

A very classical result of the theory of viscosity solution for Hamilton-Jacobi-Bellman (HJB) equations is that
solving the optimal control problem above is in fact equivalent to solving the Eikonal equation:

.. math::
  \left\lbrace \begin{matrix}
  |\nabla \phi (x) | = \frac{g(\rho(x))}{v(\rho(x))} && x \in \Omega\\
  \phi(x) = 0 && x \in \partial\Omega
  \end{matrix}\right.

We now claim that the direction field :math:`\vec{V}(t,x)` should be the unitary descending gradient of :math:`\phi`. If we suppose that there exists :math:`X^*_x` an optimal trajectory, i.e. :math:`\phi(x) = \int_0^{+\infty} \mathbb{1}_\Omega(X^*_x(t))g(\rho(X^*_x(t)))	extrm{d} t`, then we have

.. math::
  \vec{V}(t,x) = \dot{X}^*_x(0) = -\frac{\nabla \phi(x)}{|\nabla \phi(x)|}.

The Hughes model definition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Then the complete Hughes' model introduced in [Hug02]_ was the following:

.. math:: 
  \left\lbrace \begin{matrix}
  \partial_t \rho + \mathbf{div}(\vec{V}(t,x) v(t,x) \rho(t,x)) = 0 \\
  \vec{V}(t,x) = -\frac{\nabla \phi}{|\nabla \phi|} \\
  |\nabla \phi (t,x) | = \frac{g(\rho(t,x))}{v(\rho(t,x))} && x \in \Omega\\
  \phi(t,x) = 0 && x \in \partial\Omega.
  \end{matrix}\right.

.. note::
  Here, we chose to present the Hughes' model without any wall around or inside the domain for consiness' sake. Keep in mind that for a domain with walls and exits, both equations should be solved with mixed boundary condition i.e. Neumann non-crossing conditions on the walls and Dirichlet free-exit boundary conditions on the exits.

For a deep overview of the mathematical results regarding Hughes' model we defer to [AADF+23]_.

Model of Colombo-Garavello-Lécureux-Mercier
------------------

In [CGLM11]_, the authors introduced an alternative model for pedestrian flows.

This model is also featured in the present python package.

Numerical schemes
----------------

.. note::
  The presentation of the numerical schemes below is quoted from the introduction of Section 4.5.1 of [Gir25]_. It is recommanded reading it there as more detail are included.

The numerical scheme we propose consists in two coupled algorithms, each of them approximating, at one time step, one of the equations of Hughes' model (:math:`\rho` or :math:`V=\nabla \phi`, respectively) given the solution of the other equation. 

Mesh definition
^^^^^^^^^^^^^^^^

Let :math:`T>0`. Let :math:`J \in \mathbb{N}^*`. We discretize the interval :math:`[0,T]` as :math:`J + 1` time steps and we set:

.. math::
  \Delta t := \frac{T}{J}, \;\; t^j := j\Delta t.

We suppose that :math:`\Omega` is a polygonal domain. Let :math:`\Delta x > 0`.

Let :math:`\underline{\textrm{\raisebox{0pt}[1pt][0.5pt]{$\triangle$}}}, |\triangle| > 0`, we denote :math:`\Delta x := (\underline{\textrm{\raisebox{0pt}[1pt][0.5pt]{$\triangle$}}}, |\triangle|)`.
We consider a triangular mesh defined by a set of open triangles :math:`M_\Delta := (\mathcal{T}_n)_{1 \leq n \leq N}` such that

1. for all :math:`1\leq n \leq N`, the area of the triangle :math:`\mathcal{T}_n` denoted by :math:`|\mathcal{T}_n|` satisfies :math:`0 < |\mathcal{T}_n| \leq |\triangle|`.
2. for all :math:`1\leq n \leq N`, if we denote :math:`\mathcal{T}_n = A_nB_nC_n` then :math:`\max \{ |A_n B_n|, |A_n C_n|, |C_n B_n|\} \leq \underline{\textrm{\raisebox{0pt}[1pt][0.5pt]{$\triangle$}}}`.
3. for any :math:`1 \leq i,j \leq N`, we have :math:`\mathcal{T}_i \bigcap \mathcal{T}_j = \emptyset` and

.. math::
  \bar{\mathcal{T}}_i \bigcap \bar{\mathcal{T}}_j \neq \emptyset \Rightarrow
    &\textrm{ upon reparametrization of the points of the triangle }\\
     &\textrm{either } \bar{\mathcal{T}}_i \bigcap \bar{\mathcal{T}}_j = \{A_i\} = \{ A_j \},\\
    &\textrm{or } \bar{\mathcal{T}}_i \bigcap \bar{\mathcal{T}}_j = [A_iB_i] = [ A_j B_j].

4. :math:`\bigcup_n \bar{\mathcal{T}}_n = \bar\Omega`.

We use the following notations:

- the set of all the vertices is denoted by :math:`V_\Delta`.
- we say that two points of :math:`V_\Delta` are neighbours if there exists a triangle in :math:`M_\Delta` such that both points are vertices of the triangle.
- for any :math:`P \in V_\Delta`, we denote by :math:`\mathcal{V}(P)` the set of neighbours of :math:`P`.
- we denote by :math:`T_\Delta(P)` the set of all the triangles of :math:`M_\Delta` having :math:`P` as one of its vertices.

We also denote :math:`N := \bold{Card}(M_{\Delta x})` i.e.
:math:`M_{\Delta x} = (\mathcal{T}_n)_{1\leq n \leq N}`.

General overview of the scheme
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let the initial datum :math:`\rho_0` be lower semi-continuous (so that the Eikonal equation is well posed, see [Gir25]_). We discretize :math:`\rho_0` by defining, for any :math:`n \in \mathbb{N}`,

.. math::
  \rho^0_n = \inf_{x \in \mathcal{T}_n} \rho_0(x), \;\; \rho^0_\Delta(x) := \sum_{1\leq n \leq N} \mathbb{1}_{\mathcal{T}_n}(x) \rho^0_n.

Let :math:`\rho_\Delta : \bar\Omega \rightarrow \mathbb{R}` be constant on the triangles of :math:`M_{\Delta x}` i.e. :math:`\forall n \in \llbracket 1, N \rrbracket, \; \exists \rho_n \in \mathbb{R},` such that :math:`\forall x \in \mathcal{T}_n, \;\; \rho_\Delta(x) = \rho_n`.
Then, we define the lower semi-continuous envelope of :math:`\rho_\Delta` by:

.. math::
  \rho_{\Delta*}(x) := \min_{\begin{matrix}n \in \llbracket 1, N \rrbracket, \\ \textrm{ s.t. } x \in \widebar{\mathcal{T}_n}\end{matrix}} \rho_n.

Then, we will apply the following procedure iteratively for any :math:`j \in \llbracket 0 , J \rrbracket`.

1. We compute a numerical approximation :math:`\phi^j_\Delta` of the solution to the eikonal equation

.. math::
  \left\{ \begin{matrix} |\nabla \phi| = c((\rho^j_{\Delta})_*) &x \in \bar\Omega \setminus \mathcal{E}\\
  \phi(x) = 0 &x \in \mathcal{E}.\end{matrix}\right.

We use either a :math:`\mathbf{FMT}`, a :math:`\mathbf{FME}` (see [Gir25]_) or a third algorithm (detailled below) named the :math:`\mathbf{FMTC}` algorithm.

2. We compute :math:`V^j_\Delta` corresponding to the unit vector opposite to the gradient of :math:`\phi_\Delta^j` on each triangle :math:`(\mathcal{T}_n)_{1\leq n \leq N}`. In order to compute :math:`V^j_\Delta`, we denote :math:`\mathcal{T}_n = ABC`. As :math:`\phi_\Delta^j` is affine on :math:`ABC`, the gradient is constant and, for any :math:`n \in \llbracket 1, N \rrbracket`, we set:

.. math::
  V_n^j := \frac{1}{\mathcal{H}_{ABC}(\phi(A),\phi(B),\phi(C))}\times \frac{(\phi(B) - \phi(A))\vec{AB}^\bot - (\phi(C) - \phi(A))\vec{AC}^\bot}{\det(\vec{AC},\vec{AB})}.

See [Gir25]_ for more details about this formula. We define :math:`V^j_\Delta` by:

.. math::
  V^j_\Delta(x) := \sum_{1\leq n\leq N} \mathbb{1}_{\mathcal{T}_n}(x) V_n^j.

3. We compute a numerical approximation :math:`\rho^{j+1}_\Delta` of :math:`\rho(\Delta t,\cdot)` where :math:`\rho` is the solution to

.. math::
  \left\{ \begin{matrix} \rho_t + \div \left[ f(\rho) V^j_\Delta \right] = 0 &(t,x) \in (0,\Delta t) \times \Omega\\
  \rho = 0 &x \in \mathcal{E}\\
  f(\rho)(V^j_\Delta \cdot \vec{n}) = 0 &x \in \mathcal{W} \\
  \rho(0,x) = \rho^j_\Delta(x).
  \end{matrix}\right.

We use a finite volume scheme on :math:`M_\Delta` in order to approximate :math:`\rho`. We propose two different computations for the flux on the edges of the mesh, namely the **discontinuous flux** method and the **weighted flux** method. See below for details about the finite volume scheme and the different methods.


The finite volume scheme for the SCL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this paragraph, we present the explicit finite volume scheme we will use to compute :math:`(\rho^{j+1}_n)_{1\leq n \leq N}` knowing :math:`(\rho^{j}_n)_{1\leq n \leq N}` and :math:`(V^{j}_n)_{1\leq n \leq N}`.
We suppose that the maximal density :math:`\rho_{max}` equals :math:`1` i.e. we suppose that the flux function :math:`f` is Lipschitz, concave and such that
:math:`f(0) = f(1) = 0` and that :math:`\rho_0 \in [0,1]`.
We introduce a few additional notations.


- Let :math:`M \in \mathbb{N}` be the number of distinct edges in the mesh :math:`M_\Delta`. We denote by :math:`E_\Delta := (e_m)_{1\leq m \leq M}` the set of all the edges of the triangles of :math:`M_\Delta`.
- For any :math:`n \in \llbracket 1, N \rrbracket`, we denote by :math:`E_n` the set of all the edges of :math:`\mathcal{T}_n`.
- For any :math:`m \in \llbracket 1, M \rrbracket`, we denote by :math:`\mathcal E^m` the set of all triangles of :math:`M_\Delta` that admit :math:`e_m` as one of its edges. Notice that :math:`\bold{Card}(\mathcal{E}^m) \in \{1,2\}`.
- For any :math:`m \in \llbracket 1, M \rrbracket`, we denote by :math:`|e_m|` the geometrical length of the edge :math:`|e_m|`.
- For any :math:`m \in \llbracket 1, M \rrbracket`, for any :math:`\mathcal{T} \in \mathcal{E}^m`, we denote by :math:`\vec{n}_m(\mathcal{T})` the unit normal vector to :math:`e_m` that is pointing outward of :math:`\mathcal{T}`.
- For any :math:`\mathcal{T} \in M_\Delta` we denote :math:`\rho^j_{\mathcal{T}} := \rho^j_n` where :math:`\mathcal{T}_n = \mathcal{T}`. Analogously, we denote :math:`V^j_{\mathcal{T}} := V^j_n` where :math:`\mathcal{T}_n = \mathcal{T}`.

We also recall that, for any :math:`f:\mathbb{R} \rightarrow \mathbb{R}`, for any :math:`a,b \in \mathbb{R}`, the Godunov numerical flux is defined by:

.. math::
  \bold{God}_f(a,b) := \left\{ \begin{matrix} \min_{c \in [a,b]} f(c) &\textrm{if } a \leq b \\
  \max_{c \in [b,a]} f(c) &\textrm{if } b < a. \end{matrix}  \right.


The finite volume scheme corresponds to the following algorithm for any fixed :math:`j \in \llbracket 0, J \rrbracket`.

1. For any :math:`m \in \llbracket 1, M \rrbracket`, we denote

.. math::
  \mathcal{E}^m = \{  \mathcal{T} \} \textrm{ if } \bold{Card}(\mathcal{E}^m) =1,
  \mathcal{E}^m = \{  \mathcal{T}, \mathcal{T}' \} \textrm{ if } \bold{Card}(\mathcal{E}^m) =2.

We want to compute the flux :math:`f^j_m(\mathcal{T})` crossing the edge :math:`e_m` coming from :math:`\mathcal T`. We distinguish between the two following cases.

  - If :math:`\bold{Card}(\mathcal E^m) = 2`, we propose two different methods to compute :math:`f^j_m(\mathcal{T})`.
  
    - The **discontinuous flux** method: we use a dichotomy method to find :math:`k \in [0,1]` such that:
    
    .. math::
      \mathbf{God}_{V^j_{\mathcal{T} } \cdot \vec{n}_m(\mathcal{T})f}(\rho^j_{\mathcal T},k)-\mathbf{God}_{V^j_{\mathcal{T'} } \cdot \vec{n}_m(\mathcal{T})f}(k,\rho^j_{\mathcal T'}) = 0.
    
    Then we set:
    
    .. math::
      f^j_m(\mathcal{T}) := \bold{God}_{V^j_{\mathcal{T}} \cdot \vec{n}_m(\mathcal{T})f}(\rho^j_{\mathcal T},k).
    
    - The \textbf{weighted flux} method: first, if :math:`V^j{\mathcal{T}} \rho^j_{\mathcal{T}} + V^j{\mathcal{T'}} \rho^j_{\mathcal{T'}} \neq \vec{0}`, we define
    
    .. math::
      \vec{v}_m(\mathcal{T}) := \frac{V^j{\mathcal{T}} \rho^j_{\mathcal{T}} + V^j{\mathcal{T'}} \rho^j_{\mathcal{T'}}}{\left|V^j{\mathcal{T}} \rho^j_{\mathcal{T}} + V^j{\mathcal{T'}} \rho^j_{\mathcal{T'}}\right|}.
    
    If :math:`V^j{\mathcal{T}} \rho^j_{\mathcal{T}} + V^j{\mathcal{T'}} \rho^j_{\mathcal{T'}} = \vec{0}`, we set :math:`\vec{v}_m(\mathcal{T}) = \vec{0}`.
    Then we set:
    
    .. math::
      f^j_m(\mathcal{T}) :=  \bold{God}_{\vec{v}_m(\mathcal{T})\cdot \vec{n}_m(\mathcal{T}) f}(\rho^j_{\mathcal T},\rho^j_{\mathcal T'}).
    
  - Else if :math:`\bold{Card}(\mathcal E^m) = 1`, we have that :math:`e_m \in \partial \Omega`. Then, once again we distinguish between two cases.
  
    - If :math:`e_m \in \mathcal{E}`, we set:
    
    .. math::
      f^j_m(\mathcal{T}) := \bold{God}_{V^j_{\mathcal{T}} \cdot \vec{n}_m(\mathcal{T}) f}(\rho^j_{\mathcal T},0).
    
    - Else if :math:`e_m \in \mathcal{W}`, we set:
    
    .. math::
      f^j_m(\mathcal{T}) := 0.
  


2. For any :math:`n \in \llbracket 1, N \rrbracket`, we set

.. math::
  \rho^{j+1}_n := \rho^{j}_n - \frac{\Delta t}{|\mathcal{T}_n|} \sum_{e_m \in E_n} f_m^j(\mathcal{T}_n)  |e_m|.



.. note::
  The **discontinuous flux** method correspond to the algorithm detailled in \cite{AbrahamPreprint}. Heuristically, we are solving the SCL while assuming there is a discontinuity of :math:`V^j_\Delta` at each edge of :math:`E_\Delta`.
  If :math:`V_\Delta^j` is constant across the edge :math:`e_m` then the **discontinuous flux** method is equivalent to setting
  
  .. math::
    f^j_m(\mathcal{T}) := \bold{God}_{V^j_{\mathcal{T}} \cdot \vec{n}_m(\mathcal{T})f}(\rho^j_{\mathcal T},\rho^j_{\mathcal T'}),
  
  which corresponds exactly to the classical Godunov numerical flux in the continuous case. The **discontinuous flux** method correponds exactly to the use of the identity :math:`k \mapsto k` ``transmission map'', in the terminology of \cite{CancesAndreianov2015}.

.. note::
  The **weighted flux** method corresponds to the following procedure: we define :math:`\vec{v}_m(\mathcal{T})` as a weighted mean of the two vectors :math:`V^j_{\mathcal{T}}` and :math:`V^j_{\mathcal{T}'}`;
  then \eqref{eq:defMidVecFlux} corresponds to the classical Godunov numerical flux as if we had :math:`V^j_{\mathcal{T}}=V^j_{\mathcal{T}'}=\vec{v}_m(\mathcal{T})`.\\
  The heuristics for this method comes from the following situation.
  .. image:: assets/MidVectorSchema.png
  Notice that, here, we have:
  
  .. math::
    V^j_{\mathcal{T}'}\cdot \vec{n}_m(\mathcal{T}) = 0.
  
  Then, if we use the **discontinuous flux** method, for any :math:`\rho_{\mathcal{T}}^j, \rho_{\mathcal{T}'}^j  \in [0,1]` we have :math:`f_m^j(\mathcal{T}) = 0`. Then there is no flux exiting the cell :math:`\mathcal{T}` even if :math:`\rho_{\mathcal{T}'}^j = 0`.
  Heuristically, this means that the agents of cell :math:`\mathcal{T}` are prevented from moving in the direction :math:`V^j_{\mathcal{T}}` because the ``phantom'' agents of cell :math:`\mathcal{T}'` (since the cell is almost empty) should move in the incompatible direction :math:`V^j_{\mathcal{T}'}`.
  This heuristic is our inspiration to consider, as a practical alternative, the **weighted flux** method. Indeed, if we use the **weighted flux** method, if :math:`\rho_{\mathcal{T}}^j \gg \rho_{\mathcal{T}'}^j` then we have that :math:`\vec{v}_m(\mathcal{T}) \simeq V_{\mathcal{T}}^j`.
  Then
  
  .. math::
  f_m^j \simeq \max_{c \in [\rho_{\mathcal{T}'}^j,\rho_{\mathcal{T}}^j]} f(c) \;\; V_{\mathcal{T}}^j \cdot \vec{n}_m(\mathcal{T}) \gg 0.
  
  This represents a kind of ``majority-rule'' where the direction of the high density cells prevails over the direction of the low density cells.
  
.. warning::
  Even if we do not provide any proof of convergence for the above finite volume scheme in [Gir25]_, we can still derive the CFL condition that guarantees the monotonicity and the stability of the scheme. Here the CFL takes the following form:
  
  .. math::
    \Delta t \leq \frac{\underline{|\triangle|}}{3\underline{\textrm{\raisebox{0pt}[1pt][0.5pt]{$\triangle$}}}Lip_f},
  
  where :math:`\underline{|\triangle|}` denotes the minimal area of a triangle in :math:`M_\Delta`. The numerical experiments with the present package must be done under the above CFL condition !


Scheme the eikonal equation: the :math:`\mathbf{FMTC}` algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this paragraph, we fix :math:`j \in \llbracket 0, j \rrbracket`. Then, for any :math:`n \in \llbracket 0, N\rrbracket`, we set:

.. math::
  c_n := c(\rho^j_n).

We define :math:`\phi_\Delta : V_\Delta \rightarrow \mathbb{R}^+` as the solution of the :math:`\mathbf{FMTC}` numerical scheme.
Let :math:`P^0 = V_\Delta \bigcap E`. We set :math:`\phi_\Delta(P) = 0` for any :math:`P \in P^0`.
For any :math:`m \in \mathbb{N}`, we introduce the following iterative :math:`\mathbf{FMTC}` procedure, using the conceptual framework of Fast Marching algorithms where :math:`P^m` denotes the given set of ``frozen'' vertices of the mesh at step :math:`m`.

1. We compute the set of neighbours among which we will choose the next vertex (or vertices) to freeze (named the ``narrow band''):

.. math::
  NB^m = \left\{ P \in V_\Delta \setminus P^m \textrm{ s.t. } \exists Q \in P^m, P \in \mathcal{V}(Q) \right\}.

2. If :math:`NB^m = \emptyset`, the algorithm is terminated.
3. Else, for any :math:`A \in NB^m`, we compute :math:`\mathcal{V}_A` as described below.\\
  We consider :math:`(\mathcal{T}_k)_{1\leq k \leq K}` the set of all triangles of :math:`M_\Delta` such that :math:`A` is one of the vertices of the triangle and at least one other vertex of the triangle is in :math:`P^m`. Denote by :math:`ABC` the triangle :math:`\mathcal{T}_k`.\\

  For each :math:`k` we compute :math:`V^k_A` distinguishing between three cases.
    
  a. If only :math:`B` (resp. only :math:`C`) is in :math:`P_m`, then :math:`V_A^k = \phi_\Delta(B) + |AB| c_k` (resp. :math:`V_A^k = \phi_\Delta(C) + |AC| c_k` ).

  b. If :math:`B` and :math:`C` are both in :math:`P^m`, we suppose, up to renaming the vertices, that :math:`\phi_\Delta(B) \leq \phi_\Delta(C)`. We denote :math:`V_B := \phi_\Delta(B)` and :math:`V_C = \phi_\Delta(C)`. Then we set :math:`V_A^k = \tilde{V}_A` where :math:`\tilde{V}_A` is defined by \eqref{eq:defTildeVa} i.e.

  .. math::
    V_A^k := \left\{
    \begin{matrix}
        &V_B + F|AB| \qquad\qquad\qquad\qquad \textrm{ if } c_k\; \vec{AB}\cdot\vec{BC} +(V_C -V_B)|AB| >0 \\
        \\
        &V_C + F|AC| \qquad\qquad\qquad\qquad \textrm{ if } c_k\; \vec{AC}\cdot\vec{BC} +(V_C -V_B)|AC| <0 \\
        \\
        &V_B + \frac{\vec{AB}\cdot\vec{CB}}{BC^2}(V_C-V_B) + |\det(\vec{AB},\vec{CB})| \frac{\sqrt{(c_k)^2 BC^2 - (V_C-V_B)^2}}{BC^2} \qquad \textrm{ else. }\end{matrix}\right.
    
  Then we set 

  .. math::
    \mathcal{V}_A := \min_k \{ V_A^k \}.

4. We freeze the point :math:`P = \argmin_{A \in NB^m} \mathcal{V}_A` i.e. :math:`P^{m+1} = P^m \cup \{ P \}`. We set :math:`\phi_\Delta(P) = \mathcal{V}_P`. If multiple points satisfy :math:`\phi_\Delta(P) = \mathcal{V}_P`, we freeze all these points. We loop back to step 1.


In the following, we denote by :math:`\phi_\Delta` the unique function of :math:`W^{1,\infty}(\bar\Omega)` that is affine on each triangle and
:math:`\phi_\Delta(P)` is given by the above algorythm for any :math:`P \in V_\Delta`.

Bibliography
------------

.. [AADF+23] Debora Amadori, Boris Andreianov, Marco Di Francesco, Simone Fagioli, Théo Girard, Paola Goatin, Peter Markowich, Jan F. Pietschmann, Massimiliano D. Rosini, Giovanni Russo, Graziano Stivaletta, and Marie-Therese Wolfram. The mathematical theory of Hughes’ model : a survey of results. In Crowd Dynamics, Volume 4 : Analytics and Human Factors in Crowd Modeling, Modeling and Simulation in Science, Engineering and Technology. Springer, 2023.

.. [CGLM11] Rinaldo M. Colombo, Mauro Garavello, and Magali Lécureux-Mercier. A class of nonlocal models for pedestrian traffic. Mathematical Models and Methods in Applied Sciences, 22:1150023, 2011.

.. [Gir25] Théo Girard. On singular Hamilton-Jacobi-Bellman and scalar conservation laws encountered in pedestrian models. Mathematics [math]. Université de Tours, 2025.

.. [Hug02] Roger L. Hughes. A continuum theory for the flow of pedestrians. Transportation Research Part B-methodological, 36:507–535, 2002.

.. [LW55] Michael J Lighthill and Gerald Beresford Whitham. On kinematic waves. II. A theory of traffic flow on long crowded roads. Proceedings of the Royal Society of London A: Mathematical, Physical and Engineering Sciences, 229(1178):317–345, 1955.

.. [Ric56] Paul I Richards. Shock waves on the highway. Operations research, 4(1):42–51, 1956.
