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
tests can be run and the test coverage can be calculated with the following command::

    $ python -m unittest discover -v tests
    $ coverage run --source=hermipy -m unittest discover -v tests

A package automating the installation for users of the `Arch Linux`_ distribution is available on the `Arch User Repository`_,
under the name *python-hermipy-git*.
For users of other Linux distributions,
a `guix`_ package is included in the *pkgs* directory,
enabling the reproducible build of the C++ component at the core of Hermipy.
Currently, the *Boost* package provided by `guix`_ does not include a shared library file for *Boost NumPy*;
therefore, a statically-linked version of *Boost NumPy* is distributed with Hermipy in `cpp/boost`,
on which `guix`_ relies to build the C++ component.
With the `guix`_ package manager available,
buliding and testing the package can be achieved via the following commands::

    $ mkdir $HOME/hermipy && cd $HOME/hermipy
    $ wget https://raw.githubusercontent.com/urbainvaes/hermipy/master/pkgs/hermipy.scm
    $ export GUIX_PACKAGE_PATH=$HOME/hermipy
    $ guix package --install python-hermipy
    $ # Update PYTHONPATH

.. _Boost: https://en.wikipedia.org/wiki/Boost_(C%2B%2B_libraries)
.. _CMake: https://en.wikipedia.org/wiki/CMake
.. _NumPy: https://en.wikipedia.org/wiki/NumPy
.. _SciPy: https://en.wikipedia.org/wiki/SciPy
.. _SymPy: https://en.wikipedia.org/wiki/SymPy
.. _matplotlib: https://en.wikipedia.org/wiki/Matplotlib
.. _Arch Linux: https://en.wikipedia.org/wiki/Arch_Linux
.. _Arch User Repository: https://aur.archlinux.org/packages/python-hermipy-git/
.. _guix: https://www.gnu.org/software/guix/
