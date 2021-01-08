function pval = dPOLYVALsph(deg,coeff,Z,R,jvec,dbox)

%--------------------------------------------------------------------------
% Object:
% This routine evaluates the weighted-regression polynomial of degree
% "deg", defined through its coefficient vector "coeff", at the d-dim
% target points "Z".
%--------------------------------------------------------------------------
% Input:
% deg: polynomial degree
% coeff: 1-column array of of weighted-regression coefficients;
% Z: d-column array of target point coordinates;
% R: invertible triangular factor provided by dPOLYFIT;
% jvec: vector of column indices, selects a polynomial basis;
% dbox: d-dim box where to adapt the product Chebyshev basis.
%--------------------------------------------------------------------------
% Output;
% pval: evaluation of the weighted-regression polynomial at the target
%    points
%--------------------------------------------------------------------------
% Dates:
% Written on 26/07/2020 by M. Dessole, F. Marcuzzi, M. Vianello.
%
% Last update:
% 03/01/2021 by A. Sommariva.
%--------------------------------------------------------------------------
% .........................  Function Body ................................

[V,dbox] = dCHEBVAND(deg,Z,dbox);
pval = (V(:,jvec)/R)*coeff;