

function leb = dLEBsph(deg,nodes,w,X,jvec,dbox)

%--------------------------------------------------------------------------
% Object:
% This routine computes the Lebesgue constant on "nodes" of d-variate
% weighted least-squares polynomial fitting at "X".
%--------------------------------------------------------------------------
% Input:
% deg: polynomial degree;
% nodes: d-column array of points on which the Lebesgue constant is
%    required, at degree "deg";
% w: 1-column array of nonnegative weights or nonnegative scalar in
%    case of equal weights;
% X: d-column array of control point coordinates, useful to estimate the
%    Lebesgue constant;
%--------------------------------------------------------------------------
% Output:
% leb: Lebesgue constant estimate on "nodes", based on evaluations on "X".
%--------------------------------------------------------------------------
% Dates:
% Written on 26/07/2020 by M. Dessole, F. Marcuzzi, M. Vianello.
%
% Modified by:
% 29/10/2020: M. Vianello;
% 03/11/2020: M. Dessole, M. Vianello;
% 28/11/2020: A. Sommariva.
%--------------------------------------------------------------------------


% .........................  Function Body ................................

% ..... troubleshooting .....
if nargin < 6, dbox=[]; end
if nargin < 5, jvec=[]; end
if nargin < 4, X=[]; end
if nargin < 3, w=[]; end
if isempty(X), X=nodes; end
if isempty(w), w=1; end



% ........................ main code below ................................

if isempty(dbox), dbox = boxdef([X; nodes]); end

[UX,~] = dCHEBVAND(deg,X,dbox);
[VY,jvec,~,R,~] = dORTHVANDsph(deg,nodes,w,jvec,[],dbox);
VX = UX(:,jvec)/R;

% ..... Lebesgue constant approximation .....

VYVXt = VY*VX';

if isscalar(w)
    V = w*VYVXt;
else
    V = zeros(size(VYVXt));
    for k=1:size(V,2)
        V(:,k)=VYVXt(:,k).*w;
    end
end

leb = norm(V,1);



