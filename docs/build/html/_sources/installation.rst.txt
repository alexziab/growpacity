Installation
============

This page explains how to install the **growpacity** package.

Requirements
------------

- Python 3.9 or later
- NumPy
- Astropy
- ``optool`` (`Github <https://github.com/cdominik/OpTool>`_)

Optional dependencies
----------------------
- Numba (for fast interpolation)
- SciPy (as an alternative for interpolation)
- C compiler (for building the C extension)

Installation
------------

You can install the package via GitHub. If you already have ``optool`` installed system-wide, you can just clone normally:

.. code-block:: bash

    git clone https://github.com/alexziab/growpacity.git
    cd growpacity
    pip install .

Otherwise, clone with submodules to include `optool`:

.. code-block:: bash

    git clone --recurse-submodules https://github.com/alexziab/growpacity.git
    cd growpacity
    pip install .

If you already cloned without submodules but want to fetch `optool`, run the following command from within the `growpacity` directory:

.. code-block:: bash

    git submodule init && git submodule update

Installing ``optool``
---------------------

To install ``optool``, navigate to the ``optool`` subdirectory and follow the build instructions in its README. Typically, this involves running:

.. code-block:: bash

    cd growpacity/optool
    make multi=true

This will build the executable, which must be accessible to the ``growpacity`` package. This can be done by adding the ``optool`` directory to your system's ``PATH`` environment variable, or by copying the executable to a directory that is already in your ``PATH``.