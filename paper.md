---
title: "Hughes2d: approximating solutions of Hughes' model"
tags:
  - python
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

`hughes2d` is an open-source python package for simulation pedestrian crowds in two dimensions. More specifically, the package is designed to compute approximations of Hughes' model introduced in [@Hug02] (and other related models). The Hughes model is a macroscopic model -there is no agents here, the crowd is represented by a density function- coupling two non-linear partial differential equations.


# Statement of need

The mathematical modeling of pedestrian crowd is a rapidly developping topic since a few decades. There exist multiple software for simulation crowds of pedestrians both open source ([@vadere],[@JuPedSim],[@UMANS],[@Cromosim],[@CrowdWalk]) or not (PTV viswalk, MassMotion). We defer to the [awesome-crowdynamics repository](https://github.com/pozapas/awesome-crowdynamics) for an exhaustive list of the available softwares. However, up to our knowledge, all these softwares deal with microscopic simulations. We propose here a python package for macroscopic simulations of pedestrian evacuations, specifically for Hughes' model which is one of the most famous macroscopic pedestrian flow models.

# State of the field

The Hughes model has been thoroughly studied during the last two decades (see [@survey]) but there exists, at the moment, no general mathematical result of existence of solutions in 2D for this model. Some simulations appear in a few papers (see [@Goatin2014]) but in a slightly modified context. The present package should provide a reliable and open-source solution to approximate the behavior of Hughes' model.
We also hope that this package will help formulating conjectures in the future.


# Mathematical introduction to Hughes' model

The model consists in a system of two equations set on a bounded domain $\Omega \subset \mathbb{R}^2$. The first PDE is a vector-directed **scalar conservation law** with discontinuous flux. The second is an **Eikonal equation** with discontinuous source term.


> **Remark:** We defer to the introduction of Section 4.5.1 of [@Gir25] for the presentation of the numerical schemes used in the present package.

## The scalar conservation law

The first equation models the flow of a density $\rho$. To be more precise  $\rho(t,x)$ represents the density of pedestrian at time $t$ and location $x$. It is bounded between $0$ and a given $\rho_{max} >0$. We assume that the pedestrian move with the speed of agents $v(t,x)$ at time $t$ and location $x$ following a unitary direction field $\vec{V}(t,x) \in \mathcal{S}_1$.
We chose the speed of agents given by $v(\rho(t,x))$ where $v$ is a decreasing function defined on $[0,\rho_{max}]$ such that
$$ v(0)=: v_{\max} \textrm{ and } v(\rho_{max})= 0.$$
A classical example of such a speed function is $v(\rho) = v_{max}\frac{\rho_{max}-\rho}{\rho_{max}}$. This corresponds to the very classical Lighthill-Whitham-Richards (LWR) model for traffic flows (see [@LW55], [@Ric56]).

## The Eikonal equation

The second equation of the model characterizes the unitary direction field $\vec{V}(t,x)$ depending on the density $\rho$ in the whole domain. We assume that the pedestrians want to minimize their exit time while also trying to avoid high density regions. In order to model this situation, we use an optimal control problem. We suppose that the density $\rho(\cdot) \in \mathcal{C}^1(\bar \Omega)$ stays constant in time (this assumption is quite controversial).
Then, for any $x \in \Omega$, we define the value function by
$$ \phi(x) = \inf_{X \in \mathcal{A}_x} \int_0^{+\infty} \mathbb{1}_{\Omega}(X(t))g(\rho(X(t))) \textrm{d} t,$$
where $A_x$ denotes the set of admissible trajectories and $g$ is an increasing function modeling the discomfort of a agent standing in a high density region.
A very classical result of the theory of viscosity solution for Hamilton-Jacobi-Bellman (HJB) equations is that
solving the optimal control problem above is in fact equivalent to solving the Eikonal equation:


$$
 \left\lbrace \begin{matrix}
 |\nabla \phi (x) | = \frac{g(\rho(x))}{v(\rho(x))} && x \in \Omega\\
 \phi(x) = 0 && x \in \partial\Omega
 \end{matrix}\right.
 $$

## The Hughes model definition

We define the direction field $\vec{V}(t,x)$ by chosing the unitary descending gradient of $\phi$.
$$ \vec{V}(t,x) = \dot{X}^\star_x(0) = -\frac{\nabla \phi(x)}{|\nabla \phi(x)|}.$$
Then the complete Hughes' model introduced in [@Hug02] was the following:

$$
 \left\lbrace \begin{matrix}
 \partial_t \rho + \mathbf{div}(\vec{V}(t,x) v(t,x) \rho(t,x)) = 0 \\
 \vec{V}(t,x) = -\frac{\nabla \phi}{|\nabla \phi|} \\
 |\nabla \phi (t,x) | = \frac{g(\rho(t,x))}{v(\rho(t,x))} && x \in \Omega\\
 \phi(t,x) = 0 && x \in \partial\Omega.
 \end{matrix}\right.
$$

> **Remark:** Here, we chose to present the Hughes' model without any wall around or inside the domain for consiness' sake. Keep in mind that for a domain with walls and exits, both equations should be solved with mixed boundary condition i.e. Neumann non-crossing conditions on the walls and Dirichlet free-exit boundary conditions on the exits.


For a deep overview of the mathematical results regarding Hughes' model we defer to [@survey].

> **Remark:** We also included a numerical scheme for different models in the Hughes2d package such as the model of Colombo-Garavello-Lécureux-Mercier (see [@CGLM11]) but for conciseness' sake, we only present the Hughes model.

# Software design

The `hughes2d` package provides an easy way to produce simulation of macroscopic pedestrian models without having to learn the mathematical theory related to the PDEs involved in the models. Indeed a particular attention must be dedicated to the choice of numerical schemes for the PDEs involved in Hughes' model (i.e. an adapted finite volume scheme for discontinuous scalar conservation law and a fast marching algorithm for the eikonal equation). Additionnaly, most of the prototypes of code produced by mathematicians work only on square meshes. Here we chose to provide a numerical scheme working on a triangular grid in order which is much more flexible towards applications.

# Examples

We present a simple simulation of Hughes' model in a room with exits located at the end of two corridors. More precisely, we use the setting:
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
 We also mention that many videos corresponding to various numerical experiments (including this one) are available at [https://theorgirard.github.io/simulations](https://theorgirard.github.io/simulations).

![Simple Hughes simulation](docs/source/assets/demo.png)

In the above simulation, we can observe what seems to be distinctive features of Hughes' model.

- **Repartition of the agents between the different exits.** Notice that after the time $t=10s$ the agents seem to be separated in two different groups, one for each exit.

- **Geometry of the congestion figures.** Notice that after $t=15s$ the density profiles don't seem to evolve much. The room evacuates at a slow pace and the density profile for different times "look alike" until the end of the evacuation.

- **Regularization of the density in time.** Starting at $t=1s$, we can observe that the density seems to be continuous on the boundary of the support of the initial datum. This regularization in time could be a key property in order to give a mathematical proof of existence of a solution.

# Research impact statement

On a mathematical point of view, the simulation of Hughes' model suggests that there indeed exist solutions to this problem (a mathematical fact that is not proven at the moment).
Apart from the mathematical conjectures we derive from the simulations, we also hope that, by making an open source software, more researchers will take an interest in macroscopic pedestrian models. A numerical comparison of the simulation of `Hughes2d` with other models is already on-going (see [BOUM project](https://conferences.cirm-math.fr/3512.html)). 

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing
of this manuscript, or the preparation of supporting materials.


# Acknowledgements
This research was partially funded by l’Agence Nationale de la Recherche (ANR), project ANR-22-CE40-0010 COSS. I want to thank deeply Vincent Perrollaz for motivating me to publish this software and for acting as a beta tester. I also want to thank Carine Lucas for counseling me about open-source licences and standards.

# References
