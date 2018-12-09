(define-module (hermipy)
  #:use-module (guix licenses)
  #:use-module (guix packages)
  #:use-module (guix download)
  #:use-module (guix build-system python)
  #:use-module (gnu packages cmake)
  #:use-module (gnu packages python)
  #:use-module (gnu packages boost))

(define-public
  python-hermipy
  (package
    (name "python-hermipy")
    (version "v0.3.1")
    (source
      (origin
       (method url-fetch)
       (uri (string-append "https://github.com/urbainvaes/hermipy/archive/" version ".tar.gz"))
       (sha256 (base32 "1mzj5szj1fpqg66zx5ksyqqbsiwy57p1zskpsri9v5q5yxqm0lwp"))))
    (build-system python-build-system)
    (native-inputs
      `(("cmake" ,cmake)))
    (inputs
      `(("boost" ,boost)))
    (propagated-inputs
      `(("python-numpy" ,python-numpy)
        ("python-scipy" ,python-scipy)
        ("python-sympy" ,python-sympy)))
    (home-page "http://www.github.com/urbainvaes/hermipy")
    (synopsis "Python library for the Hermite Galerkin method")
    (description "Hermipy is a thin Python library that allows automating most
of the operations involved in implementing a Hermite spectral method. The
library uses the data structures provided by NumPy, and it also depends on
SymPy for symbolic calculations.")
    (license gpl3)))
