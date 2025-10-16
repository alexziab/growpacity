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
        amin=0.1, # minimum grain size [microns]
        amax_min=0.1, # smallest maximum grain size [microns]
        amax_max=1.0, # largest maximum grain size [microns]
        q_min=-4.5, # smallest power-law index
        q_max=-2.5, # largest power-law index
        Nq=5, # number of power-law indices to sample
        Namax=3, # number of maximum grain sizes to sample
        T_min=1, # lowest temperature [K]
        T_max=1000, # highest temperature [K]
        NT=100, # number of temperatures to sample
        dirc="./data", # directory to store opacity files
        optool_args="", # additional arguments to pass to OpTool (e.g., composition)
    )

We now compute and store the opacities.
In principle, you only need to run this once per set of parameters;
the results will be saved in the specified directory and can be reused later.

.. code-block:: python

    OC.execute_optool(quiet=True) # run OpTool to compute frequency-dependent opacities
    OC.build_mean_opacities() # calculate Rosseland and Planck mean opacities
    OC.compute_and_store_master_arrays() # store results in arrays

We can now access the computed opacity tables with:

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


Notes
-----

- All arguments marked with [...] can be astropy Quantities with suitable units.
- See the docstrings for each class and function for more detailed usage.
