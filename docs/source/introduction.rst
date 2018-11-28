Introduction
============

*Hermipy* is a *Python* library for the numerical solution of partial differential equations (PDEs) using the Hermite Galerking method.
The choice of the *Python* language was motivated by the following reasons:

- *Software freedom*: The reference implementation of *Python* is open source,
  with the important practical implication that
  the ability to run a *Python* program distributed under a free licence is
  not dependent upon the continual existence of a third-party organization.

- *Ecosystem*: There exist a number of mature *Python* libraries for computational mathematics,
  such as NumPy_, SciPy_ and SymPy_, which make the *Python* language particularly well-suited for the development of mathematical software.
  More recently,
  some work has also been invested in the development of the Numba_ library,
  which brings just-in-time (JIT) compilation to *Python*.

- *Extensibility*: It is possible to write *C* or *C++* extensions for *Python*,
  which can be used to improve the efficiency of a *Python* program when appropriate.
  With `Boost.Python`_ and its recent NumPy extension,
  the exchange of data between *Python* and *C++* or *C* is relatively straightforward and efficient.


When the development of *Hermipy* started,
the only *Python* library known to the author for automating spectral methods in *Python* was shenfun_,
which exposes an interface similar to that of FEniCS_,
a collection of free software components that enables automating the finite element method.
The focus of shenfun_ is on solving PDEs in bounded boxes,
by means of tensor products of Fourier, Chebyshev or Legendre bases.
In contrast, *Hermipy* aims at providing the tools for automating spectral methods in unbounded domains,
using bases of Hermite polynomials or Hermite functions.

.. _NumPy: https://en.wikipedia.org/wiki/NumPy
.. _SciPy: https://en.wikipedia.org/wiki/SciPy
.. _SymPy: https://en.wikipedia.org/wiki/SymPy
.. _Numba: http://numba.pydata.org/
.. _FEniCS: https://en.wikipedia.org/wiki/FEniCS_Project
.. _shenfun: https://github.com/spectralDNS/shenfun
.. _Boost.Python: https://www.boost.org/doc/libs/1_68_0/libs/python/doc/html/
