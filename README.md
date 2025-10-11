# Growpacity
Growpacity is a Python and C toolkit for calculating dust opacities in astrophysical environments. It provides a master class (`OpacityCalculator`) for generating, storing, and interpolating mean opacities (Rosseland and Planck) over customizable grids of grain size distributions, temperature, and power-law exponents.
This package is designed as a wrapper around [OpTool](https://github.com/cdominik/optool), allowing users to specify grain composition and obtain temperature-, maximum grain size-, and dust size distribution-dependent mean opacities. The resulting tables are lightweight and optimized for use in radiation hydrodynamics and dust coagulation models, where grain size distributions evolve dynamically.

**Note:** Growpacity does not introduce a new dust opacity model; it provides a computationally efficient framework for tabulating and interpolating opacities based on user-defined choices of grain composition and size distribution. The applicability of the results depends on the user's physical assumptions.

### Core Functionality

- **Python API:**  
    - Use `OpacityCalculator` to automate the generation of opacity tables by sweeping over grain size (`amin`, `amax`), power-law exponent (`q`), and temperature (`T`) ranges.
    - Stores absorption coefficients and mean opacities in well-defined files for later use.
    - Provides fast trilinear interpolation (`evaluate_mean_opacity`) for arbitrary parameter values.

- **C API:**  
    - High-performance routines for loading precomputed opacity arrays and evaluating mean opacities via trilinear interpolation.
    - Example usage and benchmarking included in `growpacity.c`.

### Typical Workflow

1. **Generate Opacity Tables:**  
     Instantiate `OpacityCalculator` in Python, set desired ranges, and run `execute_optool()` to generate absorption coefficients using OpTool.
2. **Compute Mean Opacities:**  
     Use `build_mean_opacities()` to calculate Rosseland and Planck means and save them.
3. **Tabulate Results:**  
     Combine results into arrays for fast lookup and interpolation.
4. **C Interpolation:**  
     Use the C routines for rapid evaluation in performance-critical applications.

### File Outputs

- `{dirc}/dustkappa_{name}_a{amax}_q{q}.inp` — Absorption coefficients
- `{dirc}/kappaRP_{name}_a{amax}_q{q}.dat` — Mean opacities
- `{dirc}/q.dat`, `amax_um.dat`, `T_K.dat` — Grid definitions
- `{dirc}/kR_cm2g.dbl`, `kP_cm2g.dbl` — Opacity arrays (binary)

### Example

See `example.ipynb` for a step-by-step demonstration of generating and visualizing opacities.

### Requirements

- Python 3, NumPy, Astropy, Numba
- [OpTool](https://github.com/cdominik/optool) (external executable)
- GCC (for compiling C code)

### Related Publication

For scientific details and methodology, see the related publication:  
**A computationally efficient dust opacity model suitable for coagulation models**  
_Alexandros Ziampras, Ludwig-Maximilians-Universität München_  
The method enables efficient tabulation and interpolation of Rosseland and Planck mean opacities for evolving dust size distributions, suitable for use in radiation hydrodynamics and coagulation models. Growpacity is not a new opacity model, but a tool for efficiently handling user-defined dust properties.  
[arXiv link or DOI if available]

### References

- For scientific details, see comments in `growpacity.py` and `growpacity.c`.
- OpTool: [https://github.com/cdominik/optool](https://github.com/cdominik/optool)
- See publication for further context and usage recommendations.

