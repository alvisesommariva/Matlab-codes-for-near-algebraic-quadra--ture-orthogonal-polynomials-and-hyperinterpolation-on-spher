
function [coeff,R,jvec,dbox] = dPOLYFITsph(deg,nodes,w,f,jvec,dbox)

%--------------------------------------------------------------------------
% Object:
% This routine computes the coefficients of the weighted-regression
% polynomial for total-degree "deg", w.r.t. a "w"-orthogonal basis on a
% d-dim point cloud "nodes".
%--------------------------------------------------------------------------
% Input:
% deg: total degree of the algebraic regression polynomial;
% nodes: d-column array of point coordinates;
% w: 1-column array of nonnegative weights, or nonnegative scalar in
%    case of equal weights (in case of doubts set it as "[]");
% f: 1-column array of sampled function values at "nodes";
% * jvec: vector of column indices that selects a polynomial basis;
% * dbox: variable that defines a hyperrectangle with sides parallel to the
%    axis, containing the domain (or pointset "nodes" in the discrete
%    case).
%    If "dbox" is not provided, it is the smaller "hyperrectangle", with
%    sides parallel to the cartesian axes, containing the pointset "nodes".
%    It is a matrix with dimension "2 x d", where "d" is the dimension of
%    the space in which it is embedded the domain.
%    For instance, for a 2-sphere, it is "d=3", for a 2 dimensional
%    polygon it is "d=2".
%    As example, the set "[-1,1] x [0,1]" is described as
%                          "dbox=[-1 0; 1 1]".
%
% Note: the variables with an asterisk "*" are not mandatory and can be
% also set as empty matrix.
% Note: it is a variant of "dPOLYFIT" suitable for the dimension of
% total degree polynomial space on the sphere.
%--------------------------------------------------------------------------
% Output:
% coeff: 1-column array of weighted regression coefficients;
% R: triangular factor in the economy size QR decomposition
%                   diag(sqrt(w))*C(:,jvec)=Q*R
%   where "C=dCHEBVAND(n,Y)" (or a particular basis if for the domain some
%   orthogonal basis are available);
% jvec: vector of column indices, selects a polynomial basis;
% dbox: variable that defines the d-dim box where to adapt the
%       product-Chebyshev basis.
%--------------------------------------------------------------------------
% Dates:
% Written on 26/07/2020 by M. Dessole, F. Marcuzzi, M. Vianello.
%
% Last update:
% 05/11/2020: A. Sommariva.
%--------------------------------------------------------------------------

if nargin < 5, jvec=[]; end
if nargin < 6, dbox=[]; end

% .............................  main code ................................

[~,jvec,Q,R,dbox]=dORTHVANDsph(deg,nodes,w,jvec,[],dbox);
coeff=Q'*(sqrt(w).*f);





