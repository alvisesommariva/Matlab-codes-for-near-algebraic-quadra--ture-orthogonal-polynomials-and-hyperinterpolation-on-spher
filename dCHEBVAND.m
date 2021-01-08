function [V,dbox] = dCHEBVAND(deg,X,dbox)

%--------------------------------------------------------------------------
% Object:
% This routine computes the Chebyshev-Vandermonde matrix for degree "deg"
% on a d-dimensional point cloud "X".
%
% The Chebyshev basis is the tensorial Chebyshev basis of total degree
% "deg", shifted on the hyperrectangle defined by "dbox".
%
% If "dbox" is not provided, the routine sets that variable to define the
% smaller "hyperrectangle" (box) with sides parallel to the cartesian
% axes and containing the pointset "X".
%--------------------------------------------------------------------------
% Input:
% deg: polynomial degree;
% X: d-column array of "m" points cloud (matrix "m x d");
% * dbox: variable that defines the smallest hyperectangle with sides
%     parallel to the axis, containing the domain.
%     If "dbox" is not provided, it defines the smaller "hyperrectangle",
%     with sides parallel to the cartesian axes, containing the pointset
%     "X".
%     It is a matrix with dimension "2 x d", where "d" is the dimension
%     of the space in which it is embedded the domain.
%     For instance, for a 2-sphere, it is "d=3", while for a 2 dimensional
%     polygon it is "d=2".
%     As example, the set "[-1,1] x [0,1]" is described as "[-1 0; 1 1]".
%
% Note: the variables with an asterisk "*" are not mandatory and can be
% also set as empty matrix.
%--------------------------------------------------------------------------
% Output:
% V: shifted Chebyshev-Vandermonde matrix for degree "deg" on the pointset
%    "X", relatively to "dbox".
% dbox: variable that defines the hyperrectangle with sides parallel to the
%     axis, containing the domain.
%--------------------------------------------------------------------------
% Data:
% The original routine has been written by M. Dessole, F. Marcuzzi and
% M. Vianello on 11/06/2020.
% It has been modified on:
% 22/10/2020 by A. Sommariva;
% 29/10/2020 by M. Vianello;
% 05/11/2020 by A. Sommariva.
%--------------------------------------------------------------------------



% ........................... Function body ...............................



% ...... troubleshooting ......

% box containing the cloud
if nargin < 3, dbox=[]; end
if isempty(dbox)
    a=min(X); b=max(X); dbox=[a;b];
else
    a=dbox(1,:); b=dbox(2,:);
end



% ..... main code below .....

% d-uples of indices with sum less or equal to "deg" graded lexicographical
% order
d=size(X,2);
N = nchoosek(deg+d,d); duples = zeros(N,d);
for i=2:N
    duples(i,:) = mono_next_grlex(d,duples(i-1,:));
end

% mapping the mesh in the hypercube "[-1,1]^d"
map = zeros(size(X));
for i=1:d
    map(:,i)=(2*X(:,i)-b(i)-a(i))/(b(i)-a(i));
end

% Chebyshev-Vandermonde matrix on the mesh
T=chebpolys(deg,map(:,1));
V=T(:,duples(:,1)+1);
for i=2:d
    T=chebpolys(deg,map(:,i));
    V=V.*T(:,duples(:,i)+1);
end










function T=chebpolys(deg,x)

%--------------------------------------------------------------------------
% Object:
% This routine computes the Chebyshev-Vandermonde matrix on the real line
% by recurrence.
%--------------------------------------------------------------------------
% Input:
% deg: maximum polynomial degree
% x: 1-column array of abscissas
%--------------------------------------------------------------------------
% Output:
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

T=zeros(length(x),deg+1);
t0=ones(length(x),1); T(:,1)=t0;
t1=x; T(:,2)=t1;

for j=2:deg
    t2=2*x.*t1-t0;
    T(:,j+1)=t2;
    t0=t1;
    t1=t2;
end









function x = mono_next_grlex ( m, x )

%*****************************************************************************80
%
%% MONO_NEXT_GRLEX: grlex next monomial.
%
%  Discussion:
%
%    Example:
%
%    M = 3
%
%    #  X(1)  X(2)  X(3)  Degree
%      +------------------------
%    1 |  0     0     0        0
%      |
%    2 |  0     0     1        1
%    3 |  0     1     0        1
%    4 |  1     0     0        1
%      |
%    5 |  0     0     2        2
%    6 |  0     1     1        2
%    7 |  0     2     0        2
%    8 |  1     0     1        2
%    9 |  1     1     0        2
%   10 |  2     0     0        2
%      |
%   11 |  0     0     3        3
%   12 |  0     1     2        3
%   13 |  0     2     1        3
%   14 |  0     3     0        3
%   15 |  1     0     2        3
%   16 |  1     1     1        3
%   17 |  1     2     0        3
%   18 |  2     0     1        3
%   19 |  2     1     0        3
%   20 |  3     0     0        3
%
%    Thanks to Stefan Klus for pointing out a discrepancy in a previous
%    version of this code, 05 February 2015.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    05 February 2015
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer X(M), the current monomial.
%    The first item is X = [ 0, 0, ..., 0, 0 ].
%
%    Output, integer X(M), the next monomial.
%

%
%  Ensure that 1 <= M.
%
if ( m < 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'MONO_NEXT_GRLEX - Fatal error!' );
    fprintf ( 1, '  M < 1\n' );
    error ( 'MONO_NEXT_GRLEX - Fatal error!' );
end
%
%  Ensure that 0 <= XC(I).
%
for i = 1 : m
    if ( x(i) < 0 )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'MONO_NEXT_GRLEX - Fatal error!' );
        fprintf ( 1, '  X(I) < 0\n' );
        error ( 'MONO_NEXT_GRLEX - Fatal error!' );
    end
end
%
%  Find I, the index of the rightmost nonzero entry of X.
%
i = 0;
for j = m : -1 : 1
    if ( 0 < x(j) )
        i = j;
        break
    end
end
%
%  set T = X(I)
%  set X(I) to zero,
%  increase X(I-1) by 1,
%  increment X(M) by T-1.
%
if ( i == 0 )
    x(m) = 1;
    return
elseif ( i == 1 )
    t = x(1) + 1;
    im1 = m;
elseif ( 1 < i )
    t = x(i);
    im1 = i - 1;
end

x(i) = 0;
x(im1) = x(im1) + 1;
x(m) = x(m) + t - 1;

return



