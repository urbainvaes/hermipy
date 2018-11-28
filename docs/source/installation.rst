Installation
============

Hermipy comes with a C++ component that needs to be compiled before the module can be used.
To compile it, the following dependencies are required:

- The Boost_ library, including the *Boost Python* and *Boost Numpy* components;
- The CMake_ cross-platform build system.

With the dependencies present,
the C++ component can be compiled by running the following command from the root directory of the repository::

    $ python setup.py build

The module can then be installed::

    $ python setup.py install --prefix ~/.local

To install globally, prepend ``sudo`` to this command and omit the prefix.
In addition to the two dependencies above,
at runtime the module depends on the following python packages:

- NumPy_ and SciPy_, for efficient scientific computing routines;
- SymPy_, for symbolic calculations;
- matplotlib_, for plots and visualization.

Once installed,
tests can be run with the following command::

    $ python -m unittest discover -v tests

A package automating the installation for users of the `Arch Linux`_ distribution is available on the `Arch User Repository`_,
under the name *python-hermipy-git*.

.. _Boost: https://en.wikipedia.org/wiki/Boost_(C%2B%2B_libraries)
.. _CMake: https://en.wikipedia.org/wiki/CMake
.. _NumPy: https://en.wikipedia.org/wiki/NumPy
.. _SciPy: https://en.wikipedia.org/wiki/SciPy
.. _SymPy: https://en.wikipedia.org/wiki/SymPy
.. _matplotlib: https://en.wikipedia.org/wiki/Matplotlib
.. _Arch Linux: https://en.wikipedia.org/wiki/Arch_Linux
.. _Arch User Repository: https://aur.archlinux.org/packages/python-hermipy-git/
