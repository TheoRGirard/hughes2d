---
title: "hughes2d: Approximating solutions to the Hughes pedestrian-flow model"
tags:
  - Python
  - non-linear discontinuous PDE
  - pedestrian evacuation
  - macroscopic crowd dynamics
  - finite volume scheme
  - fast marching method
authors:
  - name: Théo Girard
    orcid: 0000-0002-9382-6746
    affiliation: 1
affiliations:
  - name: Institut Denis Poisson, Université de Tours, France
    index: 1
date: 02 december 2025
bibliography: paper.bib
---

# Summary

`hughes2d` is an open-source python package for simulation pedestrian crowds in two dimensions. More specifically, the package is designed to compute approximations of the Hughes model introduced in [@Hug02] (and other related models). The Hughes model is a macroscopic model -- there are no particles here, the crowd is represented by a density function- coupling two non-linear partial differential equations.


# State of the field

The mathematical modeling of pedestrian crowd is a rapidly developing topic since a few decades. There exist multiple software packages for simulating crowds of pedestrians -- both open source [@vadere;@JuPedSim;@UMANS;@Cromosim;@CrowdWalk]) and proprietary (PTV viswalk, MassMotion). We defer to the [awesome-crowdynamics repository](https://github.com/pozapas/awesome-crowdynamics) for an exhaustive list of the available packages. However, to the best of our knowledge, all these developments deal with microscopic simulations. We propose here a Python package for macroscopic simulations of pedestrian evacuations, specifically for the seminal Hughes pedestrian flow model.


# Statement of need

The Hughes model has been thoroughly studied during the last two decades [see @survey], but there exists, at the moment, no general mathematical result of existence of solutions in 2D for this model. While example numerical solutions to the Hughes model can be found in literature (e.g., [@Goatin2014]), these are not supported by a publicly available and reusable codebase. The present package is aimed at providing a reliable and open-source solution to approximate the behavior of the Hughes model.
We also hope that this package will help formulating conjectures in the future.


# Software design

The `hughes2d` package provides a way to produce simulation of macroscopic pedestrian models without having to learn the mathematical theory related to the underlying PDEs. Indeed, for the Hughes model, particular attention must be dedicated to the choice of numerical schemes (i.e., an adapted finite volume scheme for discontinuous scalar conservation law and a fast marching algorithm for the eikonal equation). Additionally, instead of the common choice of square meshes, we extend the range of applications by employing a numerical scheme working on a triangular grid.


# Mathematical introduction to the Hughes model

The model consists of a system of two equations set on a bounded domain $\Omega \subset \mathbb{R}^2$. The first PDE is a vector-directed **scalar conservation law** with discontinuous flux. The second is an **eikonal equation** with a discontinuous source term.


> **Remark:** We defer to the introduction of Section 4.5.1 in @Gir25 for the presentation of the numerical schemes used in the present package.

## The scalar conservation law

The first equation models the flow of a density $\rho(t,x)$ representing the density of pedestrians at time $t$ and location $x$. It is bounded between $0$ and a given $\rho_{\max} >0$. We assume that the pedestrian move with the speed of agents $v(t,x)$ at time $t$ and location $x$ following a unitary direction field $\vec{V}(t,x) \in \mathcal{S}_1$.
We chose the speed of agents given by $v(\rho(t,x))$ where $v$ is a decreasing function defined on $[0,\rho_{\max}]$ such that
$$ v(0)=: v_{\max} \textrm{ and } v(\rho_{\max})= 0.$$
A classical example of such a speed function is $v(\rho) = v_{\max}\frac{\rho_{\max}-\rho}{\rho_{\max}}$. This corresponds to the very classical Lighthill-Whitham-Richards (LWR) model for traffic flows [see @LW55;@Ric56].

## The eikonal equation

The second equation of the model characterizes the unitary direction field $\vec{V}(t,x)$ depending on the density $\rho$ in the whole domain. We assume that the pedestrians want to minimize their exit time while also trying to avoid high density regions. In order to model this situation, we use an optimal control problem. We suppose that the density $\rho(\cdot) \in \mathcal{C}^1(\bar \Omega)$ stays constant in time (this assumption is quite controversial).
Then, for any $x \in \Omega$, we define the value function by
$$ \phi(x) = \inf_{X \in \mathcal{A}_x} \int_0^{+\infty} \mathbb{1}_{\Omega}(X(t))g(\rho(X(t))) \textrm{d} t,$$
where $A_x$ denotes the set of admissible trajectories and $g$ is an increasing function modeling the discomfort of a agent standing in a high density region.
A very classical result of the theory of viscosity solution for the Hamilton-Jacobi-Bellman (HJB) equations is that
solving the optimal control problem above is in fact equivalent to solving the eikonal equation:


$$
 \left\lbrace \begin{matrix}
 |\nabla \phi (x) | = \frac{g(\rho(x))}{v(\rho(x))} && x \in \Omega\\
 \phi(x) = 0 && x \in \partial\Omega
 \end{matrix}\right.
 $$

## The Hughes model definition

We define the direction field $\vec{V}(t,x)$ by choosing the unitary descending gradient of $\phi$.
$$ \vec{V}(t,x) = \dot{X}^\star_x(0) = -\frac{\nabla \phi(x)}{|\nabla \phi(x)|}.$$
Then the complete model introduced in @Hug02 is defined by:

$$
 \left\lbrace \begin{matrix}
 \partial_t \rho + \mathbf{div}(\vec{V}(t,x) v(t,x) \rho(t,x)) = 0 \\
 \vec{V}(t,x) = -\frac{\nabla \phi}{|\nabla \phi|} \\
 |\nabla \phi (t,x) | = \frac{g(\rho(t,x))}{v(\rho(t,x))} && x \in \Omega\\
 \phi(t,x) = 0 && x \in \partial\Omega.
 \end{matrix}\right.
$$

> **Remark:** Here, we chose to present the Hughes model without any wall around or inside the domain for the sake of brevity. Keep in mind that for a domain with walls and exits, both equations should be solved with a mixed boundary condition, i.e., Neumann non-crossing conditions on the walls and Dirichlet free-exit boundary conditions on the exits.


For an overview of the mathematical results regarding the Hughes model, we refer the reader to @survey.

> **Remark:** the `hughes2d` package also includes a numerical scheme for other models, including the Colombo-Garavello-Lécureux-Mercier model [see @CGLM11].


# Examples

We present a simple simulation using the Hughes model in a room with exits located at the end of two corridors. More precisely, we use the setting:
$$
    \bar\Omega := [-2,0]\times [3,4] \cup [0,10]\times [0,7] \cup [10,12] \times [3,4],
    $$
    $$
    \mathcal{E} := \lbrace -2\rbrace \times [3,4] \cup \lbrace 12\rbrace \times [3,4], \;\;
  \mathcal{W} := \partial \Omega \setminus \mathcal{E},
  $$
  $$
    \rho_0(x) := 0.7\times \mathbb{1}_{B((7,2.5),2.4)}, \; \;
    c(\rho) = 1 +5\rho.
$$
Videos corresponding to various numerical experiments (including this one) are available at [https://theorgirard.github.io/simulations](https://theorgirard.github.io/simulations).

![Sample simulation using the Hughes model](docs/source/assets/demo.png)

In the above simulation, the following distinctive features of the Hughes model can be observed:

- **Repartition of the agents between the different exits.** Notice that after the time $t=10$s, the agents seem to be separated in two different groups, one for each exit.

- **Geometry of the congestion figures.** Notice that after $t=15$s, the density profiles do not seem to evolve much. The room evacuates at a slow pace and the density profile for different times "look alike" until the end of the evacuation.

- **Regularization of the density in time.** Starting at $t=1$s, we can observe that the density seems to be continuous on the boundary of the support of the initial datum. This regularization in time could be a key property in order to give a mathematical proof of existence of a solution.

# Research impact statement

On a mathematical point of view, the simulation using the Hughes model suggests that there indeed exist solutions to the problem (a mathematical fact that is not proven at the moment).
Apart from the mathematical conjectures we derive from the simulations, we also hope that, by making an open source software, more researchers will take an interest in macroscopic pedestrian models. A numerical comparison of the simulation using `hughes2d` with other models is already on-going (see [BOUM project](https://conferences.cirm-math.fr/3512.html)).

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing
of this manuscript, or the preparation of supporting materials.


# Acknowledgments
This research was partially funded by l’Agence Nationale de la Recherche (ANR), project ANR-22-CE40-0010 COSS. I want to thank deeply Vincent Perrollaz for motivating me to publish this software and for acting as a beta tester. I also want to thank Carine Lucas for counseling me about open-source licenses and standards.

# References
