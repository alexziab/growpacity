# Growpacity
`growpacity` is a Python and C toolkit for calculating dust opacities in astrophysical environments. It provides a framework for generating, storing, and interpolating mean opacities (Rosseland and Planck) over customizable grids of grain size distributions, temperature, and powerlaw exponents.

This package is designed as a wrapper around [optool](https://github.com/cdominik/optool), allowing users to specify grain composition and obtain temperature-, maximum grain size-, and dust size distribution-dependent mean opacities. The resulting tables are lightweight and optimized for use in radiation hydrodynamics and dust coagulation models, where grain size distributions evolve dynamically.

**Note:** Growpacity does not introduce a new dust opacity model; it provides a computationally efficient framework for tabulating and interpolating opacities based on user-defined choices of grain composition and size distribution. The applicability of the results depends on the user's physical assumptions.

<!-- table of contents with hyperlinks: -->
- [Installation](#installation)
     - [1. Downloading the module](#1-downloading-the-module)
     - [2. Building `optool`](#2-building-optool)
     - [3. Installing `growpacity`](#3-installing-growpacity)
- [Core Functionality](#core-functionality)
- [Typical Workflow](#typical-workflow)
- [File Outputs](#file-outputs)
- [Example](#example)
- [Requirements](#requirements)
     - [Mandatory](#mandatory)
     - [Optional but recommended](#optional-but-recommended)
- [Related Publication](#related-publication)
- [References](#references)

### Core Functionality

- **Python API:**  
    - Uses `OpacityCalculator` to generate opacity tables by sweeping over grain size (`amin`, `amax`), power-law exponent (`q`), and temperature (`T`) ranges.
    - Stores absorption coefficients and Rosseland/Planck mean opacities in ASCII or binary files for later use.
    - Provides fast trilinear interpolation (`evaluate_mean_opacity`) for arbitrary parameter values.

- **C API:**  
    - High-performance routines for loading precomputed opacity arrays and evaluating mean opacities via trilinear interpolation.
    - Example usage and benchmarking included in `growpacity.c`.

### Installation

#### 1. Downloading the module

`growpacity` ships with `optool` as a git submodule. To install, clone the repository with submodules:

```bash
git clone --recurse-submodules https://github.com/alexziab/growpacity.git
cd growpacity
```

_note: if you already cloned without submodules, run:_
```bash
git submodule init && git submodule update
```
to fetch `optool`.

#### 2. Building `optool`

After that, install `optool` by following the instructions in the `optool/` directory. Typically, this will look like:

```bash
cd growpacity/optool
make multi=true
```

Important: Ensure that the `optool` binary is in your system's PATH for `growpacity` to function correctly (i.e., you can run the command `optool` from any directory). This can be done by running `make install bindir=~/bin/` from the `optool` directory, adjusting the path as needed, or adding it to your PATH manually:

```bash
export PATH=$PATH:/path/to/optool
```

#### 3. Installing `growpacity`

Finally, install the `growpacity` Python package by navigating to the parent directory and running:

```bash
pip install .
```

Done correctly, the file structure should look like this:

```growpacity/
├── growpacity/
│   ├── __init__.py
│   ├── growpacity.c
│   ├── optool/
│   │   └── ...
│   └── ...
├── pyproject.toml
├── README.md
├── docs/
│   └── ...
├── tests/
│   └── ...
└── ...
```

### Typical Workflow

1. **Generate Opacity Tables:**  
     Instantiate `OpacityCalculator` in Python, set desired ranges, and run `execute_optool()` to generate absorption coefficients using optool.
2. **Compute Mean Opacities:**  
     Use `build_mean_opacities()` to calculate Rosseland and Planck means and save them.
3. **Tabulate Results:**  
     Combine results into arrays for fast lookup and interpolation.
4. **C Interpolation:**  
     Use the C routines for fast evaluation in performance-critical applications such as hydrodynamics simulations.

### File Outputs

- Absorption coefficients ($\kappa_\text{abs}$, $\kappa_\text{sca}$, $g$)
- Mean opacities ($\kappa_\text{R}$, $\kappa_\text{P}$) for each parameter combination
- `q.dat`, `amax_um.dat`, `T_K.dat` — Grid definitions
- `kR_cm2g.dbl`, `kP_cm2g.dbl` — Opacity arrays (binary)

### Example

See `tests/example.ipynb` for a step-by-step demonstration of generating and visualizing opacities.

_note: the example notebook requires `matplotlib` for plotting, and (optionally) `scipy` to compare the `numba` implementation against._

### Requirements

#### Mandatory

- Python `>=3.9`, Numpy, Astropy
- [optool](https://github.com/cdominik/optool), installed and accessible via terminal (see [installation-2](#2-building-optool))

#### Optional but recommended

- Numba (for JIT acceleration of Python interpolation)
- Scipy (for interpolation if Numba is not used, e.g. with python >=3.14 as of Nov 2025)
- GCC (for compiling the C extension)

### Related Publication

For scientific details and methodology, see the related publication:  
**A computationally efficient dust opacity model suitable for coagulation models**  
_Alexandros Ziampras, Tilman Birnstiel_
[arXiv link to be added]

### References

- Sphinx documentation can be found in the `docs/` folder.
- For additional details, see docstrings/comments in `growpacity/growpacity.py` and `growpacity/growpacity.c`.
- optool: [https://github.com/cdominik/optool](https://github.com/cdominik/optool)
- See publication for further context and usage recommendations.