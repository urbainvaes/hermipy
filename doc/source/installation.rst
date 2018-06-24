Installation
============

Hermipy comes with a C++ component that needs to be compiled before the module can be used.
To compile it, the following dependencies are required:

- The Boost_ library, including the *Boost python* and *Boost numpy* components;
- The CMake_ cross-platform build system.

.. _Boost: https://en.wikipedia.org/wiki/Boost_(C%2B%2B_libraries)
.. _CMake: https://en.wikipedia.org/wiki/CMake

With the dependencies present,
the C++ component can be compiled by running the following command from the root directory of the repository::

    $ python setup.py build

The module can then be installed::

    $ python setup.py install --prefix ~/.local

To install globally, prepend ``sudo`` to this command and omit the prefix.

In addition to the two dependencies above,
the module depends on the following python packages:

- NumPy_ and SciPy_, for efficient scientific computing routines;
- SymPy_, for symbolic calculations;
- matplotlib_, for plots and visualization.

.. _NumPy: https://en.wikipedia.org/wiki/NumPy
.. _SciPy: https://en.wikipedia.org/wiki/SciPy
.. _SymPy: https://en.wikipedia.org/wiki/SymPy
.. _matplotlib: https://en.wikipedia.org/wiki/Matplotlib
