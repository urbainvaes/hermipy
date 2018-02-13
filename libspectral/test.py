import spectral
def function(v): return 1
quadrature = spectral.Quad([100,100,100])
print(quadrature.integrate(function))
