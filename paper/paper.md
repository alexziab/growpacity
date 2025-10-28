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
    affiliation: [1,2]
  - name: Tilman Birnstiel
    orcid: 0000-0002-1899-8783
    affiliation: [1,3]
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

`growpacity` is a Python and C toolkit for calculating dust opacities in astrophysical environments. It provides a framework for generating, storing, and interpolating Rosseland and Planck mean opacities over customizable grids of grain size distributions, temperatures, and powerlaw exponents.
This package is designed as a wrapper around [OpTool](https://github.com/cdominik/optool), allowing users to specify a grain composition of their choice and obtain tabulated mean opacities. The resulting tables are lightweight and optimized for use in radiation hydrodynamics and dust coagulation models, where grain size distributions evolve dynamically.

# Statement of need

Radiation hydrodynamics simulations often require a dust opacity model under the assumption of a particular grain size or size distribution. It is however not straightforward to implement a dust opacity model that accounts for the dynamical evolution of this grain size distribution. This shortcoming is especially relevant in environments where dust coagulation is important, such as protoplanetary disks.

Dust coagulation models evolve the dust grain size distribution as a function of position and time, with grain sizes ranging from $\sim0.1$\,\mu$m to $\sim$cm—12 orders of magnitude in mass. The complexity of modeling dust dynamically in this manner is exacerbated by the vastly different opacities between small and large grains, with small grains dominating the opacity budget but also being subject to rapid depletion due to coagulation. A dust opacity model that can account for these effects is therefore essential in accurately capturing the thermal structure of dusty astrophysical environments such as protoplanetary disks, where the temperature profile directly influences processes such as planet formation, disk chemistry, and observational signatures of substructure. With recent strides in dust coagulation modeling [@StammlerBirnstiel:2022; e.g., @PfeilEtal:2024; @RobinsonEtal:2024], there is a growing need for an accurate, computationally efficient, flexible, and easy-to-implement dust opacity model that can be integrated into such simulations.

With `growpacity`, we outline a method to construct such a model and provide a `python` package that interfaces to `OpTool` to compute temperature-, maximum grain size-, and dust size distribution-dependent mean opacities, tabulated in a lightweight format and complemented with efficient interpolation methods for usage in radiation hydrodynamics simulations.

# Method

For a given grain composition and assuming a grain size distribution with number density $n(a)\propto a^q$ for grains with size $a$ between minimum and maximum grain sizes $a_{\min}$ and $a_{\max}$,
the [OpTool](https://ui.adsabs.harvard.edu/abs/2021ascl.soft04010D) package [@DominikEtal:2021] can compute the absorption and scattering opacities $\kappa_\text{abs}(\nu)$ and $\kappa_\text{sca}(\nu)$ (in cm$^2$/g) as well as the asymmetry factor $g(\nu)$ [@HenyeyGreenstein:1941] over a frequency grid $\nu$. This calculation is done using the Distribution of Hollow Spheres method [DHS, @MinEtal:2005]. The Rosseland and Planck mean opacities $\kappa_\text{R}$ and $\kappa_\text{P}$ can then be computed as

$$
\kappa_\text{P}(T) = \frac{\int_0^\infty \kappa_\text{abs} B_\nu(T) d\nu}{\int_0^\infty B_\nu(T) d\nu}, \quad \kappa_\text{R}(T) = \frac{\int_0^\infty u_\nu(T) d\nu}{\int_0^\infty [\kappa_\text{abs}+(1-g)\,\kappa_\text{sca}]^{-1} u_\nu(T) d\nu}, \quad u_\nu(T) = \left.\frac{dB_\nu}{dT}\right|_T,
$$

where $B_\nu(T)$ is the Planck function at temperature $T$. By fixing the grain composition and $a_{\min}$, `OpTool` can be used to calculate the absorption and scattering opacities over a grid of $q$ and $a_{\max}$. The above equation can then be used to compute and tabulate $\kappa_\text{R}$ and $\kappa_\text{P}$ over $q$, $a_{\max}$, and $T$. The resulting 3D tables can be used in any context where $q$ and $a_{\max}$ are dynamically evolved according to a dust coagulation model [@BirnstielEtal:2017; @StammlerBirnstiel:2022; @PfeilEtal:2024; @RobinsonEtal:2024].

# Interpolation algorithm

Once the 3D tables of $\kappa_\text{R}(q,a_{\max},T)$ and $\kappa_\text{P}(q,a_{\max},T)$ have been computed, they can be interpolated within the range of tabulated values. We interpolate for $\log\kappa$ as a function of $q$, $\log a_{\max}$, and $\log T$, with regular sampling in this space (i.e., logarithmic spacing for $a_{\max}$ and $T$). This works best when the mean opacities follow a powerlaw with respect to temperature $\kappa\propto T^b\Rightarrow \log\kappa\propto b\log T$, which is a reasonable approximation and especially holds for small grains [@BellLin:1994; @SemenovEtal:2003]. It also ensures that the interpolated opacities are always positive in case extrapolation is needed.

We use a trilinear interpolation scheme, which is fast and simple to implement, and takes advantage of the fact that the arrays $q$, $\log a_{\max}$, and $\log T$ are sorted and regularly spaced to efficiently locate the indices of the grid points that surround the point of interest. For each array $x\in\{q,\log a_{\max},\log T\}$ and for a target value $x_t$, we first find the index $i$ such that $x_i\leq x_t < x_{i+1}$ as $i = \lfloor(x_t-x_0)\,\Delta x^{-1}\rfloor$, where $x_0$ is the first (smallest) value in the sampling space and $\Delta x = x_{i+1}-x_{i}$ is the (constant) sampling spacing. As this information is known \emph{a priori}, this reduces the complexity of finding the required indices from $\mathcal{O}(\log N)$ to $\mathcal{O}(1)$, where $N$ is the number of grid points in $x$, and is especially efficient for larger grids. Of course, care must be taken to ensure that the indices are within the bounds of the grid.

# Implementation

The above method is implemented in the `growpacity` package. The package provides a simple `python` interface to `OpTool` to compute accurate dust opacities for a given grain composition (provided in `OpTool` format) and for a fixed $a_{\min}$, over a grid of $a_{\max}$, $q$, and $T$. The resulting 3D tables of $\kappa_\text{R}$ and $\kappa_\text{P}$ are saved in binary files that can be easily loaded with `python`. We also provide an implementation of the interpolating function in `C`, that can be readily used in radiation hydrodynamics codes. These files are especially lightweight: a typical calculation involving $q\in\left[-4.5,-2.5\right]$ with $\Delta q = 0.25$, $a_{\max}\in\left[0.1\,\mu\text{m},1\,\text{m}\right]$ sampled twice per decade, and $T\in\left[1, 2000\right]\,\text{K}$ with 100 points results in two arrays of $9\times15\times100$ elements, or about 220 kB of memory. This choice of spacing means that opacities are exactly evaluated for $q\in[-3.75,-3.5,-3]$, corresponding to dust size distributions in equilibrium due to small/large-scale turbulence, radial drift, or in a non-equilibrium growth regime [for more information, see @PfeilEtal:2024].

The code and documentation are available at  
**<https://github.com/alexziab/growpacity>**.

# Limitations and extensions

We underscore that our intent is not to provide a new or more realistic dust opacity model, but rather one suitable for use in coagulation models, where dust densities and distributions can vary as a function of position and time---the applicability of the model depends on entirely user-defined choices. The method can be easily extended to include gas opacities [e.g., @MalyginEtal:2014; @SemenovEtal:2003] and prescriptions for the sublimation of dust species [e.g., @IsellaNatta:2005], for a more complete opacity model in regimes where gas opacities are significant.

# Acknowledgements

AZ and TB acknowledge funding from the European Union under the Horizon Europe Research and Innovation Programme 101124282 (EARLYBIRD). Views and opinions expressed are those of the authors only.

# References