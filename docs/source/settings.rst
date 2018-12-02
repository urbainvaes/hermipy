Configuration
=============

A dictionary of options can be modified to configure *Hermipy*.
The default options that apply when the module is loaded are the following::

    settings = {
            'cache': False,
            'cachedir': '.cache',
            'debug': False,
            'sparse': False,
            'tensorize': True,
            'trails': False,
            }

Changing the value of an option is achieved by modifying the appropriate dicionary entry.
For example, the code below loads *Hermipy*,
enables the caching feature,
and sets the cache directory to */tmp/hermipy-cache/*::

    import hermipy as hm
    hm.settings['cache'] = True
    hm.settings['cachedir'] = '/tmp/hermipy-cache'

The available settings are described below.

``cache`` and ``cachedir``
      Setting ``cache`` to ``True`` enables the caching feature of *Hermipy*;
      with this option, the results of all the low-level calls to the C++ compiled component are saved for later reuse.
      These include in particular the results of Hermite transforms and of the computation of stiffness matrices.
      This feature works by hashing the arguments of function calls to test for availability of the result in the cache,
      and by using the *NumPy* functions ``save`` and ``load`` to write to and read from the cache.
      It is most useful when used in combination with the ``tensorize`` option.

``tensorize``
    In many applications, the transform associated to functions in several dimensions can be obtained from lower dimensional transforms by tensorization,
    reducing the computational cost significantly.
    Likewise, it is often possible to take advantage of the fact that many operators are tensorizable to speed up the computations of the corresponding stiffness matrices.
    By setting ``tensorize`` to ``True``, *Hermipy* will automatically split functions and operators to make the computations more efficient.

``sparse``
    Setting ``sparse`` to ``True`` instructs *Hermipy* to use sparse matrices for the representation of stiffness matrices.
    This is useful only when the operator being discretized admits a sparse representation in space of basis functions being considered,
    i.e. for operators with polynomial (resp. trigonometric) coefficients when Hermite polynomials (resp. Fourier series) are used.

``debug``
    Setting ``debug`` to ``True`` enables the display of debug information when using *Hermipy*.
    This is useful mainly for debugging purposes by contributors.

``trails``
    Setting ``trails`` to ``True`` enables the display of information related to the function calls performed by *Hermipy*.
    A potential use of this feature is to identify the computational bottleneck of a simulation code,
    or to monitor the progress of computationally demandinm functions.
