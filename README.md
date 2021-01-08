# MatlabCodesForNearAlgebraicQuadrature-OrthogonalPolynomials-HyperinterpolationOnSphericalTriangles

» Quadrature rules and hyperinterpolation on spherical triangles (2021) 
Matlab codes for determining quadrature rules and hyperinterpolation on spherical triangles.
The rules have (numerically) fixed algebraic degree of precision "n". Tchakaloff formulas, with a cardinality at most "(n+1)^2" can be extracted by means of a new implementation of Lawson-Hanson algorithm.
Polynomial hyperinterpolation at degree "n" is determined by means of the rules above of degree "2n" and discrete orthogonal polynomials.

» Object: Cubature rules on spherical triangles. 

The main routine is
cub_sphtri.m that determines full and compressed quadrature rules with fixed (numerical) algebraic degree of precision on spherical triangles; it requires dCATCHsph.m, dCHEBVAND.m, gqcircsect.m, LHDM.m, trigauss.m. 

The software includes the demos files:
demo_01.m that tests the cubature rules on spherical harmonics (it requires sphharmVAND.m and plot_s2.m),
demo_02.m that tests the cubature rules on a battery of tests functions. 

» Object: Polynomial hyperinterpolation on spherical triangles. 

The software includes the demos files:
demo_hyperinterpolation_01.m that tests the hyperinterpolants of several functions on some spherical triangles,
demo_hyperinterpolation_02.m that shows how to determine the hyperinterpolant on a single example.
demo_hyperinterpolation_03.m that shows estimates of the Lesbegue constant based on Christoffel function. 

» Sources: 
Papers: 
Near-algebraic Tchakaloff-like quadrature on spherical triangles (submitted)
Numerical hyperinterpolation over spherical triangles (submitted)
Last version: [MATLAB CODES (zip file)]
Older version: [MATLAB CODES (zip file)]
» Last update: January 07, 2021. 
