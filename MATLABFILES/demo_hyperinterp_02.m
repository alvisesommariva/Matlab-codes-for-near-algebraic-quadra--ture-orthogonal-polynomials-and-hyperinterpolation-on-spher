
function demo_hyperinterp_02

%--------------------------------------------------------------------------
% Object:
% Demo on algebraic polynomial hyperinterpolation on spherical triangles.
% It is a simple version to understand how to use the subroutines.
%--------------------------------------------------------------------------
% Reference paper:
% A. Sommariva and M. Vianello
% Numerical hyperinterpolation overspherical triangles
%--------------------------------------------------------------------------
% Dates:
% Written on 07/01/2021: A. Sommariva and M. Vianello.
%--------------------------------------------------------------------------

clear;

% ....................... spherical triangle vertices  ....................

% The vertices are described in cartesian coordinates as rows of "vertices"
a=0.9; b=0.5;
vertices=[0 0 1; 0 sqrt(a) sqrt(1-a); sqrt(b) sqrt(a/2) sqrt(1-a/2-b)];


% ....................... function to approximate  ........................
x0=0; y0=0; z0=1;
g=@(x,y,z) exp(-((x-x0).^2 + (y-y0).^2 + (z-z0).^2));

% ....................... hyperinterpolation degree  ......................
deg=5;

% ........................ Main code below ................................

P1=vertices(1,:); P2=vertices(2,:); P3=vertices(3,:);

[XYZW,XYZWC,momerr]=cub_sphtri(2*deg,P1,P2,P3);

% ... test pointset ...
X_test=XYZW(:,1:3); w_test=XYZW(:,4);

% ... hyperinterpolation points ...
X=XYZWC(:,1:3); w=XYZWC(:,4);

% ... evaluate function to approximate ...
g_X=feval(g,X(:,1),X(:,2),X(:,3));

% ... determine polynomial hyperinterpolant ...
[coeff,R,jvec,dbox] = dPOLYFITsph(deg,X,w,g_X);

% ... evaluate hyperinterpolant at test pointset ...
p_X_test=dPOLYVALsph(deg,coeff,X_test,R,jvec,dbox);

% ... rough estimate of hyperinterpolant error ...
g_X_test=feval(g,X_test(:,1),X_test(:,2),X_test(:,3));

RE_L2err=sqrt(w_test'*(p_X_test-g_X_test).^2)/sqrt(w_test'*(g_X_test.^2));

% ... statistics ...
fprintf('\n \t deg: %2.0f, RE_L2err: %1.2e \n \n',deg,RE_L2err);




