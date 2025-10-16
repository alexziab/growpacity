Usage
=====

This page explains how to use the **growpacity** package.

Importing the package
---------------------

After `Installation <installation.html>`_, you can simply import the package in a python environment with:

.. code-block:: python

    import growpacity

Example: building an opacity model
----------------------------------

First, we create an instance of the `OpacityCalculator` class, specifying
the desired parameters for the grain size distribution and temperature range.

.. code-block:: python

    OC = growpacity.OpacityCalculator(
        amax_max=1.0, # largest maximum grain size [microns]
        Namax=3, # number of maximum grain sizes to sample
        dirc="./data", # directory to store opacity files
    )

We now compute and store the opacities.
In principle, you only need to run this once per set of parameters;
the results will be saved in the specified directory and can be reused later.

.. code-block:: python

    OC.execute_optool(quiet=True) # run OpTool to compute frequency-dependent opacities
    OC.build_mean_opacities() # calculate Rosseland and Planck mean opacities
    OC.compute_and_store_master_arrays() # store results in arrays

The computations may take a few minutes depending on the parameters.

You will now find the following files in the specified directory:

- ``dustkappa...inp`` (ASCII):
  OpTool output files for each combination of power-law index ``q``
  and maximum grain size ``amax``. Compatible with `RADMC-3D <https://github.com/dullemond/radmc3d-2.0>`_.

- ``kappaRP...dat`` (binary):
  Contain the Rosseland and Planck mean opacities as a function of temperature,
  for each combination of power-law index ``q`` and maximum grain size ``amax``.

- ``{q, amax_um, T_K}.dat`` (ASCII):
  Text files containing the sampled values of power-law index, maximum grain size (in microns),  
  and temperature (in Kelvin). The first line in each file indicates the number of sampled values.

- ``{kR, kP}_cm2g.dbl`` (binary):
  Files containing the 3D arrays of Rosseland and Planck mean opacities,
  with shape ``(Nq, Namax, NT)``.

Example: loading and using pre-computed opacities
-------------------------------------------------

After having computed and stored the opacities once, you can load them in a new python session
without needing to recompute them, as long as you specify the same parameters to the `OpacityCalculator` class.
We can access the computed opacity tables with:

.. code-block:: python

    q      = OC.q # array of power-law indices with Nq elements
    amax   = OC.amax.to_value('micron') # array of maximum grain sizes with Namax elements
    T      = OC.T.to_value('K') # array of temperatures with NT elements
    kR, kP = OC.load_master_arrays() #each has shape (Nq, Namax, NT)

The 3D arrays ``kR`` and ``kP`` contain the Rosseland and Planck mean opacities, respectively.
You can now interpolate to get the Rosseland and Planck mean opacities for any combination of
power-law index, maximum grain size, and temperature within the specified ranges.

Note that the interpolating function expects all arguments in CGS units (i.e., amax in cm).
We still use amax in microns otherwise for convenience when interfacing to OpTool.

.. code-block:: python

    q_test = -3.65
    amax_test = 0.55 # microns
    T_test = 150 # K

    # Rosseland mean opacity
    kR_test = growpacity.evaluate_mean_opacity((q, amax, T), kR, q_test, amax_test/1e4, T_test)

    # Planck mean opacity
    kP_test = growpacity.evaluate_mean_opacity((q, amax, T), kP, q_test, amax_test/1e4, T_test)

The module also provides a convenience function that wraps ``evaluate_mean_opacity``
for array-like inputs:

.. code-block:: python

    q_tests = -3.5
    amax_tests = np.array([0.1, 0.5]) # microns
    T_tests = [100, 200, 1000] # K

    # Rosseland mean opacities
    kR_tests = growpacity.evaluate_mean_opacity_array((q, amax, T), kR, \
                            q_tests, amax_tests/1e4, T_tests)

Notes
-----

- All arguments marked with [...] can be astropy Quantities with suitable units.
- See the docstrings for each class and function for more detailed usage.
