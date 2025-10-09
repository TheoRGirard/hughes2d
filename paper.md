---
title: "Hughes2d: approximating solutions of the Hughes'model"
tags:
  - python
  - non-linear discontinuous PDE 
  - pedestrian evacuation
  - macroscopic crowd dynamics
  - finite volume scheme
  - fast marching method
authors:
  - name: Théo R. Girard
    orcid: 0000-0002-9382-6746
    equal-contrib: true
    affiliation: 1
affiliations:
  - name: Institut Denis Poisson, Université de Tours
    index: 1
    ror: 05djhd259
date: 09 October 2025
bibliography: paper.bib
---

# Summary

`hughes2d` is an open-source python package for simulation pedestrian crowds in two dimensions. More specifically, the package is designed to compute approximations of Hughes' model introduced in [Hug02]. The Hughes model is a macroscopic model -there is no agents here, the crowd is represented by a density function- coupling two non-linear partial differential equations. 


# Statement of need

The mathematical modeling of pedestrian crowd is a rapidly developping topic since a few decades. There exist multiple software for simulation crowds of pedestrians both open source ([vadere],[JuPedSim],[UMANS],[Cromosim]) or not. However, up to our knowledge, all these softwares deal with microscopic simulations. We propose here a python package for macroscopic simulations of pedestrian evacuations, specifically for Hughes' model which is one of the most famous macroscopic pedestrian flow models.

The Hughes model has been thoroughly studied during the last two decades (see [Survey]) but there exists, at the moment, no general mathematical result of existence of solutions in 2D for this model. Some simulations appear in a few papers (see [Goatin]) but in a slightly modified context.
We hope that this package will help formulating conjectures in the future.


# Short introduction to Hughes' model

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements


# References

