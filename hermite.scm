(define-module (python-hermite)
  #:use-module (gnu packages)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix build-system gnu)
  #:use-module (guix packages)
  #:use-module (guix download)
  #:use-module (guix git-download)
  #:use-module (guix utils)
  #:use-module (guix build utils)
  #:use-module (guix build-system cmake))

(define-public
  python-hermite
  (package
    (name "python-hermite")
    (version "v0.1")
    (source
      (origin
        (method url-fetch)
        (uri (string-append
               "https://github.com/urbainvaes/hermite/archive/tag/"
               version ".tar.gz"))
        (sha256
          (base32 "10y3frwich8263zj344czsvfq31fwcva758cx33g2v43y6yss9xb"))))
    (build-system python-build-system)
    (inputs
      `(("boost" ,boost)))
    (native-inputs
      `(("cmake" ,cmake)))
    (propagated-inputs
      `(("python-numpy" ,python-numpy)
        ("python-scipy" ,python-scipy)
        ("python-sympy" ,python-scipy)))
    (arguments
      '(
        #:phases
        (modify-phases
          %standard-phases
          (add-after
            'unpack 'pull-mshmet
            (lambda* (#:key source inputs #:allow-other-keys)
                     (begin
                       (zero? (system* "mkdir" "-p" "download/pkg"))
                       (zero? (system* "cp" (assoc-ref inputs "mshmet") "download/pkg/mshmet.2012.04.25.tgz"))
                       (zero? (system* "cp" (assoc-ref inputs "tetgen") "download/pkg/tetgen1.5.1-beta1.tar.gz"))
                       (zero? (system* "cp" (assoc-ref inputs "hpddm") "download/pkg/hpddm.zip"))
                       )))
          (add-before
            'build 'build-modules
            (lambda _ (begin
                        (zero? (system* "mkdir" "download/include" "download/lib"))
                        (zero? (system* "make" "-C" "download/hpddm"))
                        (zero? (system* "make" "-C" "download/mshmet"))
                        (zero? (system* "make" "-C" "download/tetgen")))))
          ; (add-before
          ;   'configure 'fix-hardcoded-paths
          ;   (lambda _ (begin
          ;               (substitute* "configure" (("/bin/sh") (which "sh"))))))
        )))
    (home-page "https://github.com/urbainvaes/hermite")
    (synopsis "Library for the Hermite spectral method")
    (description "This library provides")
    (license (list
               license:lgpl3 ; mshmet
               license:lgpl2.1+))))
