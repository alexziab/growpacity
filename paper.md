---
title: "growpacity: A computationally efficient dust opacity model suitable for coagulation models"
tags:
  - Python
  - astrophysics
  - protoplanetary disks
  - dust opacities
  - radiative transfer
authors:
  - name: Alexandros Ziampras
    orcid: 0000-0003-2966-0419
    affiliation: 1
  - name: Tilman Birnstiel
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Ludwig-Maximilians-Universität München, Universitäts-Sternwarte, Scheinerstr. 1, 81679 München, Germany
    index: 1
  - name: Max Planck Institute for Astronomy, Königstuhl 17, 69117 Heidelberg, Germany
    index: 2
  - name: Exzellenzcluster ORIGINS, Boltzmannstr. 2, 85748 Garching, Germany
    index: 3
date: 2025-10-09
bibliography: refs.bib
---

# Summary

`growpacity` provides a computationally efficient method and implementation for calculating Rosseland and Planck mean dust opacities for evolving grain size distributions.  
Dust opacities are essential in radiation hydrodynamics simulations, but incorporating dynamically changing grain size distributions is typically cumbersome.  
This package offers a fast and simple way to construct opacity tables that depend on temperature, the maximum grain size, and the dust size distribution slope, suitable for use in dust coagulation and evolution models.

# Statement of need

Radiation hydrodynamics and dust evolution models in protoplanetary disks require accurate mean opacities that depend on the local dust properties.  
Standard tabulated opacities assume fixed grain size distributions, making them unsuitable for models where the dust evolves dynamically.  
`growpacity` addresses this by coupling existing dust opacity tools to an efficient tabulation and interpolation scheme, providing a practical solution for simulations where the dust size distribution varies in space and time.

# Method

For a given grain composition and a grain size distribution with number density \( n(a) \propto a^q \) for grain sizes between \( a_{\min} \) and \( a_{\max} \), absorption and scattering opacities \( \kappa_{\text{abs}}(\nu) \) and \( \kappa_{\text{sca}}(\nu) \), as well as the asymmetry factor \( g(\nu) \), are computed using the [OpTool](https://ui.adsabs.harvard.edu/abs/2021ascl.soft04010D) package [@dominik-etal-2021].  
This calculation uses the Distribution of Hollow Spheres method [@min-etal-2005], and provides frequency-dependent opacities that are then integrated to obtain the Planck and Rosseland mean opacities:

\[
\kappa_P(T) = \frac{\int_0^\infty \kappa_{\text{abs}} B_\nu(T)\, d\nu}{\int_0^\infty B_\nu(T)\, d\nu}, \qquad
\kappa_R(T) = \frac{\int_0^\infty u_\nu(T)\, d\nu}{\int_0^\infty [\kappa_{\text{abs}} + (1-g)\kappa_{\text{sca}}]^{-1} u_\nu(T)\, d\nu},
\]
where \( u_\nu(T) = (dB_\nu/dT)_T \) and \( B_\nu(T) \) is the Planck function.

Fixing the grain composition and \( a_{\min} \), `growpacity` computes mean opacities over a 3D grid of \( a_{\max} \), \( q \), and \( T \), tabulating \( \kappa_R(q, a_{\max}, T) \) and \( \kappa_P(q, a_{\max}, T) \).  
These tables are then suitable for use in dust evolution models where \( a_{\max} \) and \( q \) evolve dynamically [@birnstiel-etal-2017; @stammler-birnstiel-2022; @pfeil-etal-2024; @robinson-etal-2024].

# Interpolation and performance

Interpolation is performed on \(\log \kappa\) in the space of \((q, \log a_{\max}, \log T)\) using a trilinear interpolation scheme.  
The sampling is uniform in this space, taking advantage of the approximately power-law dependence \(\kappa \propto T^b\), which is accurate for small grains [@bell-lin-1994; @semenov-etal-2003].  
This design ensures positive opacities even during limited extrapolation.  
Because the grid spacing is regular, the indices of bounding grid points can be computed analytically, reducing lookup complexity from \(\mathcal{O}(\log N)\) to \(\mathcal{O}(1)\).  

A typical dataset with \( q \in [-4.5, -2.5] \) (Δq = 0.25), \( a_{\max} \in [0.1\,\mu\text{m}, 1\,\text{m}] \) (two samples per decade), and \( T \in [1, 2000]\,\text{K} \) (100 samples) results in two arrays of \( 9\times15\times100 \) elements—about 220 kB total.  

# Implementation

The package is written in Python and interfaces directly with `OpTool` to compute frequency-dependent opacities, then integrates and tabulates mean opacities.  
A C implementation of the interpolation routine is included for integration into radiation hydrodynamics or coagulation codes.  
The code and documentation are available at  
**<https://github.com/alexziab/growpacity>**.

# Limitations and extensions

`growpacity` does not introduce a new opacity model; it provides a framework suitable for dynamically evolving dust size distributions.  
The method can be extended to include gas opacities [@semenov-etal-2003; @malygin-etal-2014] and dust sublimation prescriptions [@isella-natta-2005], enabling broader applicability in thermochemical disk models.

# Acknowledgements

AZ acknowledges funding from the European Union under the Horizon Europe Research and Innovation Programme 101124282 (EARLYBIRD).  
Views and opinions expressed are those of the author only.

# References
