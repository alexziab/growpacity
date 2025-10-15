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

.. code-block:: python

    import growpacity

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

    OC.execute_optool(quiet=True) # run OpTool to compute frequency-dependent opacities
    OC.build_mean_opacities() # calculate Rosseland and Planck mean opacities
    OC.compute_and_store_master_arrays() # store results in arrays

Notes
-----

- All arguments marked with [...] can be astropy Quantities with suitable units.
- See the docstrings for each class and function for more detailed usage.
