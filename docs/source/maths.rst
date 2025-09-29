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

.. latex::

  In \cite{Hughes2002ACT}, Hughes proposed a mathematical model for the two-dimensional dynamics of pedestrian crowds. The model consists in a system of two equations set on a bounded domain $\Omega \subset \mathbb{R}^2$. The first equation models the flow of a density $\rho$. To be more precise  $\rho(t,x)$ represents the density of pedestrian at time $t$ and location $x$. It is bounded between $0$ and a given $\rho_{max} >0$. We assume that the pedestrian move with the speed of agents $v(t,x)$ at time $t$ and location $x$ following a unitary direction field $\vec{V}(t,x) \in \mathcal{S}_1$.
  Then, if we write the conservation of the mass on pedestrian on each subdomain of $\Omega$, we end up with the following scalar conservation law:
  \begin{equation}\label{eq:SCLHughes}
  	\partial_t \rho + \div(\vec{V}(t,x) v(t,x) \rho(t,x)) = 0.
  \end{equation}
  In Hughes' model, we assume that the speed of agents $v(t,x)$ at time $t$ and location $x$ only depends on the density $\rho(t,x)$ and that this speed is decreasing with respect to the density.
  Then, in the following we will denote the speed by $v(\rho(t,x))$ where $v$ is a decreasing function defined on $[0,\rho_{max}]$ such that
  $$v(0)=: v_{\max} \textrm{ and } v(\rho_{max})= 0.$$
  A classical example of such a speed function is $v(\rho) = v_{max}\frac{\rho_{max}-\rho}{\rho_{max}}$. This corresponds to the very classical Lighthill-Whitham-Richards (LWR in the following) model for traffic flows (see \cite{lighthill1955kinematic,richards1956shock}).\\
  
  The second equation of the model characterizes the unitary direction field $\vec{V}(t,x)$ depending on the density $\rho$ in the whole domain. We assume that the pedestrians want to minimize their exit time while also trying to avoid high density regions. In order to model this situation, we use an optimal control problem. We suppose that the density $\rho(\cdot) \in \mathcal{C}^1(\bar \Omega)$ stays constant in time (we refer to Remark \ref{rem:eikoConstantTime} for a discussion on this modeling assumption).
  Let $x \in \Omega$. For any $\alpha(\cdot) \in \mathcal{C}^1((0,+\infty),\mathcal{S}_1)$, we say that $X^\alpha_x(\cdot)$ is a trajectory controlled by $\alpha$ starting at $x$ if $X$ is a solution to the Cauchy problem:
  \begin{equation}\label{eq:HugheseikonalDynamic}
  	\left\lbrace \begin{matrix}
  	\dot{X}^\alpha_x(s) = \alpha(s)v(\rho(t,X^\alpha_x(s))) \\
  	X^\alpha_x(0) = x
  	\end{matrix}\right.
  \end{equation}
  For any $x \in \Omega$, we denote by $\mathcal{A}_x = \{ X^\alpha_x, \alpha \in \mathcal{C}^1((0,+\infty),\mathcal{S}_1) \}$ the set of all controlled trajectories starting at $x$. We define $\phi(x)$ as the minimal exit time starting at location $x$, that is to say :
  \begin{equation}\label{eq:ExitTimeHughes}
  \phi(x) = \inf_{X \in \mathcal{A}_x} \int_0^{+\infty} \mathbb{1}_{\Omega}(X(t)) \d t.
  \end{equation}
  In Hughes' model, we also take into account the discomfort caused by being surrounded by a high density crowd. In order to model this discomfort, we introduce an increasing function $g(\rho)$ with respect to the density $\rho$. The function $g(\rho)$ can be interpreted as a running cost we are paying along a trajectory $X$ for being in high density regions. Then, equation \eqref{eq:ExitTimeHughes} becomes:
  \begin{equation}\label{eq:ValueFunctionHughes}
  \phi(x) = \inf_{X \in \mathcal{A}_x} \int_0^{+\infty} \mathbb{1}_{\Omega}(X(t))g(\rho(X(t))) \d t.
  \end{equation}
  A very classical result of the theory of viscosity solution for Hamilton-Jacobi-Bellman equations is that
  solving this optimal control problem \eqref{eq:HugheseikonalDynamic}-\eqref{eq:ValueFunctionHughes} is in fact equivalent to solving the Eikonal equation \eqref{eq:EikonalHughes}:
  \begin{equation}\label{eq:EikonalHughes}
  \left\lbrace \begin{matrix}
  |\nabla \phi (x) | = \frac{g(\rho(x))}{v(\rho(x))} && x \in \Omega\\
  \phi(x) = 0 && x \in \partial\Omega
  \end{matrix}\right.
  \end{equation}
  \begin{remark}
  Here, we chose to present the Hughes' model without any wall around or inside the domain. The walls are not relevant for the one-dimensional version of the problem. We discuss the complexity introduced by the walls in Chapter \ref{chap:Hughes2D}.
  \end{remark}
  We now claim that the direction field $\vec{V}(t,x)$ should be the unitary descending gradient of $\phi$. If we suppose that there exists $X^*_x$ an optimal trajectory, i.e. $\phi(x) = \int_0^{+\infty} \mathbb{1}_\Omega(X^*_x(t))g(\rho(X^*_x(t)))\d t$, then we have
  $$\vec{V}(t,x) = \dot{X}^*_x(0) = -\frac{\nabla \phi(x)}{|\nabla \phi(x)|}.$$
  Then the complete Hughes' model introduced in \cite{Hughes2002ACT} was the following:
  \begin{equation}\label{eq:Hughes2DsansMur}
  \left\lbrace \begin{matrix}
  \partial_t \rho + \div(\vec{V}(t,x) v(t,x) \rho(t,x)) = 0 \\
  \vec{V}(t,x) = -\frac{\nabla \phi}{|\nabla \phi|} \\
  |\nabla \phi (t,x) | = \frac{g(\rho(t,x))}{v(\rho(t,x))} && x \in \Omega\\
  \phi(t,x) = 0 && x \in \partial\Omega.
  \end{matrix}\right.
  \end{equation}

Colombo-Garavello-
------------------

Bibliography
------------

.. [AADF+23] Debora Amadori, Boris Andreianov, Marco Di Francesco, Simone Fagioli, Théo Girard, Paola Goatin, Peter Markowich, Jan F. Pietschmann, Massimiliano D. Rosini, Giovanni Russo, Graziano Stivaletta, and Marie-Therese Wolfram. The mathematical theory of Hughes’ model : a survey of results. In Crowd Dynamics, Volume 4 : Analytics and Human Factors in Crowd Modeling, Modeling and Simulation in Science, Engineering and Technology. Springer, 2023.

.. [CGLM11] Rinaldo M. Colombo, Mauro Garavello, and Magali Lécureux-Mercier. A class of nonlocal models for pedestrian traffic. Mathematical Models and Methods in Applied Sciences, 22:1150023, 2011.

.. [Hug02] Roger L. Hughes. A continuum theory for the flow of pedestrians. Transportation Research Part B-methodological, 36:507–535, 2002.
