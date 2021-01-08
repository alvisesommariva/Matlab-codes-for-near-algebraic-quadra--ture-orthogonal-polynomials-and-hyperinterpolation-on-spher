
function tw=trigauss(n,alpha,beta,method)

%--------------------------------------------------------------------------
% AUTHORS.
%--------------------------------------------------------------------------
% Alvise Sommariva and Marco Vianello, University of Padova
% July 25, 2016.
%
% previous versions with the help of G. Da Fies.
%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% It computes the angles and weights of a trigonometric quadrature formula
% on [alpha,beta], 0<beta-alpha<=pi, matching the moments up to 10^(-14).
%
% Depending on the choice of the variable 'method', some methods are
% implemented.
%
% The formula integrates the canonical trigonometric basis with accuracy
% from about 10^(-15) (small omega) to about 10^(-13) (omega-->pi)
% up to n=300.
%--------------------------------------------------------------------------
% INPUT.
%--------------------------------------------------------------------------
% n: trigonometric degree of exactness.
%
% [alpha,beta]: angular interval, 0<(beta-alpha)/2<=pi.
%
% method:
%         1. 'classic' implements classical trigauss rule.
%         2. 'legendre' implements (shifted) Gauss-Legendre rule
%             matching the trigonometric moments up to 10^(-14). Few nodes
%             if angles are small otherwise the cardinality might be higher
%             than in 'classic'.
%         3. 'better': chooses 'classic' or 'legendre', so to have the
%             minimal number of points, still matching the moments up to
%             10^(-14).
%         4.  'antigauss': antigaussian formula based on Legendre rule of  
%             degree M=M(n). 
%             It provides a (M+1) x 4 matrix whose entries are 
%                             tw=[tAGL wAGL tGLe wGLe]
%             where 
%             a) (tAGL,wAGL) is the antigaussian rule;
%             b) (tGLe(1:end-1),wGLe(1:end-1)) is the gaussian rule;
%         5.  'kronrod': Gauss-Kronrod formula based on Legendre rule of  
%             degree M=M(n). 
%             It provides a (2*M+1) x 3 matrix whose entries are 
%                             tw=[tGK wGK wGLf]
%             where 
%             a) (tGK,wGK) is the Gauss-Kronrod rule;
%             b) (tGK(2:2:end),wGLf(2:2:end)) is the gaussian rule;
%         6.  'antitrigauss': antigaussian formula based on trigonometric  
%             gaussian rule of degree n (whose cardinality is n+1). 
%             It provides a (n+2) x 4 matrix whose entries are 
%                             tw=[tAGL wAGL tGLe wGLe];
%             where 
%             a) (tAGL,wAGL) is the antigaussian rule;
%             b) (tGLe(1:end-1),wGLe(1:end-1)) is the gaussian rule;
%         7.  'trig_kronrod': Gauss-Kronrod based on trigonometric gaussian 
%             rule of degree n. 
%             It provides a (2*(n+1)+1) x 3 matrix whose entries are 
%                             tw=[tTK wTK wTf];
%             where 
%             a) (tTK,wTK) is the Gauss-Kronrod rule;
%             b) (tTK(2:2:end),wTf(2:2:end)) is the trig. gaussian rule;
%         otherwise: if no string is given it chooses 'better' option by 
%             default.
%--------------------------------------------------------------------------
% OUTPUT.
%--------------------------------------------------------------------------
% tw: 1) for 'classic', 'legendre', 'better', it is a 
%     M x 2 array of (angles,weights) (M depends on the choosen rule)
%     2) for 'antigauss' and 'antitrigauss' it is a M x 4 matrix (see
%     explanation above at points 4. and 6.)
%     3) for 'kronrod' and 'trig_kronrod' it is a M x 3 matrix (see
%     explanation above at points 5. and 7.)
%--------------------------------------------------------------------------
% EXAMPLES.
%--------------------------------------------------------------------------
%
% >> format long e;
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'better');
% >> twB
% 
% twB =
% 
%     -1.508419779244729e-02     1.590094129717467e-03
%     -1.251400776401954e-02     3.493153120682099e-03
%     -8.255043791082401e-03     4.927692470361337e-03
%     -2.881384626391015e-03     5.697023547188071e-03
%      2.881384626391019e-03     5.697023547188069e-03
%      8.255043791082396e-03     4.927692470361324e-03
%      1.251400776401955e-02     3.493153120682099e-03
%      1.508419779244729e-02     1.590094129717463e-03
% 
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'classic');
% >> twB
% 
% twB =
% 
%     -1.536597391335750e-02     8.744544324142629e-04
%     -1.393392018525980e-02     1.972635825710960e-03
%     -1.146915301262857e-02     2.926255469814080e-03
%     -8.153889677130131e-03     3.662992813352760e-03
%     -4.233938891704753e-03     4.128095270647346e-03
%      2.178667008161081e-18     4.287058912019128e-03
%      4.233938891704759e-03     4.128095270647326e-03
%      8.153889677130132e-03     3.662992813352761e-03
%      1.146915301262857e-02     2.926255469814094e-03
%      1.393392018525980e-02     1.972635825710951e-03
%      1.536597391335750e-02     8.744544324142609e-04
% 
% >> format short e
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'antigauss');
% >> twB
% 
% twB =
% 
%   -1.5612e-02   5.4060e-04  -1.5084e-02   1.5901e-03
%   -1.4036e-02   2.5826e-03  -1.2514e-02   3.4932e-03
%   -1.0564e-02   4.2828e-03  -8.2550e-03   4.9277e-03
%   -5.6645e-03   5.4042e-03  -2.8814e-03   5.6970e-03
%   -7.0636e-19   5.7956e-03   2.8814e-03   5.6970e-03
%    5.6645e-03   5.4042e-03   8.2550e-03   4.9277e-03
%    1.0564e-02   4.2828e-03   1.2514e-02   3.4932e-03
%    1.4036e-02   2.5826e-03   1.5084e-02   1.5901e-03
%    1.5612e-02   5.4060e-04            0            0
% 
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'kronrod');
% >> twB
% 
% twB =
% 
%   -1.5604e-02   2.7995e-04            0
%   -1.5084e-02   7.7659e-04   1.5901e-03
%   -1.4045e-02   1.2956e-03            0
%   -1.2514e-02   1.7537e-03   3.4932e-03
%   -1.0561e-02   2.1404e-03            0
%   -8.2550e-03   2.4607e-03   4.9277e-03
%   -5.6659e-03   2.7029e-03            0
%   -2.8814e-03   2.8494e-03   5.6970e-03
%    3.4551e-18   2.8973e-03            0
%    2.8814e-03   2.8494e-03   5.6970e-03
%    5.6659e-03   2.7029e-03            0
%    8.2550e-03   2.4607e-03   4.9277e-03
%    1.0561e-02   2.1404e-03            0
%    1.2514e-02   1.7537e-03   3.4932e-03
%    1.4045e-02   1.2956e-03            0
%    1.5084e-02   7.7659e-04   1.5901e-03
%    1.5604e-02   2.7995e-04            0
% 
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'antitrigauss');
% >> twB
% 
% twB =
% 
%   -1.5366e-02   8.7445e-04  -1.5298e-02   1.0473e-03
%   -1.3934e-02   1.9726e-03  -1.3588e-02   2.3476e-03
%   -1.1469e-02   2.9263e-03  -1.0672e-02   3.4414e-03
%   -8.1539e-03   3.6630e-03  -6.8077e-03   4.2296e-03
%   -4.2339e-03   4.1281e-03  -2.3385e-03   4.6420e-03
%    2.1787e-18   4.2871e-03   2.3385e-03   4.6420e-03
%    4.2339e-03   4.1281e-03   6.8077e-03   4.2296e-03
%    8.1539e-03   3.6630e-03   1.0672e-02   3.4414e-03
%    1.1469e-02   2.9263e-03   1.3588e-02   2.3476e-03
%    1.3934e-02   1.9726e-03   1.5298e-02   1.0473e-03
%    1.5366e-02   8.7445e-04            0            0
% 
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'trig_kronrod');
% >> twB
% 
% twB =
% 
%   -1.5640e-02   1.8370e-04            0
%   -1.5298e-02   5.1143e-04   1.0473e-03
%   -1.4611e-02   8.6012e-04            0
%   -1.3588e-02   1.1787e-03   2.3476e-03
%   -1.2265e-02   1.4628e-03            0
%   -1.0672e-02   1.7183e-03   3.4414e-03
%   -8.8397e-03   1.9398e-03            0
%   -6.8077e-03   2.1160e-03   4.2296e-03
%   -4.6243e-03   2.2427e-03            0
%   -2.3385e-03   2.3207e-03   4.6420e-03
%   -4.7437e-18   2.3475e-03            0
%    2.3385e-03   2.3207e-03   4.6420e-03
%    4.6243e-03   2.2427e-03            0
%    6.8077e-03   2.1160e-03   4.2296e-03
%    8.8397e-03   1.9398e-03            0
%    1.0672e-02   1.7183e-03   3.4414e-03
%    1.2265e-02   1.4628e-03            0
%    1.3588e-02   1.1787e-03   2.3476e-03
%    1.4611e-02   8.6012e-04            0
%    1.5298e-02   5.1143e-04   1.0473e-03
%    1.5640e-02   1.8370e-04            0
%
%--------------------------------------------------------------------------
% NOTE.
%--------------------------------------------------------------------------
% For examples about the usage of formulas of antigaussian or Kronrod type,
% see the file "demo_trigauss_error.m".
%--------------------------------------------------------------------------
%% Copyright (C) 2016
%% Gaspare Da Fies, Alvise Sommariva, Marco Vianello.
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%%
%% Authors: 
%% Gaspare Da Fies, Alvise Sommariva, Marco Vianello.
%%
%% Date: JULY 25, 2016
%--------------------------------------------------------------------------

if nargin <= 3
    method = 'classic';
end

validStrings = {'classic','legendre','better','antigauss',...
    'kronrod','antitrigauss','trig_kronrod'};

if ( any(strcmpi(method, validStrings)) )
    
    if ( strcmpi(method, 'classic') )
        tw=trigauss_classic(n,alpha,beta);
    end
    
    if ( strcmpi(method, 'legendre') )
        if beta-alpha==2*pi
            tw=circle_trigquad(n,alpha,beta);
        else
            NN=degree_trigGL(n,alpha,beta);
            tw=gauss_legendre(NN,alpha,beta);
        end
    end
    
    if ( strcmpi(method, 'better') )
        if beta-alpha==2*pi
            tw=circle_trigquad(n,alpha,beta);
        else
            NN=degree_trigGL(n,alpha,beta);
            if NN <= n
                tw=gauss_legendre(NN,alpha,beta);
            else
                tw=trigauss_classic(n,alpha,beta);
            end
        end
    end
    
    if ( strcmpi(method, 'antigauss') )
        NN=degree_trigGL(n,alpha,beta);
        tw=antigauss_full(NN,alpha,beta);
    end

    if ( strcmpi(method, 'kronrod') )
        NN=degree_trigGL(n,alpha,beta);
        tw=gauss_kronrod_full(NN,alpha,beta);
    end
    
    if ( strcmpi(method, 'antitrigauss') )
        tw=antitrigauss_full(n,alpha,beta);
    end
    
    if ( strcmpi(method, 'trig_kronrod') )
        tw=trigauss_kronrod_full(n,alpha,beta);
    end
    
else
    warning('wrong string as method, choosen classical trigauss');
    tw=trigauss(n,alpha,beta,'classic');
end





function NN=degree_trigGL(n,alpha,beta)

% this routine computes the number of nodes that Gauss-Legendre needs so to
% match the trigonometric moments up to 10^(-14).

theta=(beta-alpha)/2;

if n < 500
    
    u=[1 (25:25:500)];
    
    L=[8 29    45    60    75    89   103   119   132   146   158   173  ...
        185   199 212   225   240   253   264   279   291];
    
    s=spline(u,L);
    
    NN=ceil(ppval(s,n*theta));
    
else
    
    u=n*theta;
    s=0.54*u+21;
    NN=ceil(s);
    
end



function tw=gauss_kronrod_full(n,alpha,beta)
n0=ceil(3*n/2)+1;
ab0=r_jacobi(n0,0,0);

% KRONROD.
xw=kronrod(n,ab0);
xGK=xw(:,1); wGK=xw(:,2);
tGK=(beta+alpha)/2+(beta-alpha)*xGK/2;
wGK=(beta-alpha)*wGK/2;


% GAUSS.
xwGL=gauss(n,ab0);
xGL=xwGL(:,1); wGL=xwGL(:,2);
tGL=(beta+alpha)/2+(beta-alpha)*xGL/2;

wGLl=(beta-alpha)*wGL/2; 
wGLf=zeros(size(wGK));

wGLf(2:2:end)=wGLl;

tw=[tGK wGK wGLf];





function tw=antigauss_full(n,alpha,beta)

ab=r_jacobi(n+1,0,0);
ab(end,2)=2*ab(end,2);
xwAGL=gauss(n+1,ab);

xAGL=xwAGL(:,1); wAGL=xwAGL(:,2);
tAGL=(beta+alpha)/2+(beta-alpha)*xAGL/2;
wAGL=(beta-alpha)*wAGL/2;

xwGL=gauss(n,ab(1:end-1,:));
xGL=xwGL(:,1); wGL=xwGL(:,2);
tGL=(beta+alpha)/2+(beta-alpha)*xGL/2; tGLe=[tGL; 0];
wGL=(beta-alpha)*wGL/2; wGLe=[wGL; 0];

tw=[tAGL wAGL tGLe wGLe];







function tw=trigauss_kronrod_full(n,alpha,beta)

omega=(beta-alpha)/2;
n0=ceil(3*n/2)+1;
ab0=r_trigauss(n0,omega);

% KRONROD.
xw=kronrod(n,ab0);
xTK=xw(:,1); wTK=xw(:,2);
tTK(:,1)=2*asin(sin(omega/2)*xTK(:,1))+(beta+alpha)/2;


% GAUSS.
xwT=gauss(n,ab0);
xT=xwT(:,1); wT=xwT(:,2);
tT(:,1)=2*asin(sin(omega/2)*xT(:,1))+(beta+alpha)/2;

wTf=zeros(size(wTK));

wTf(2:2:end)=wT;

tw=[tTK wTK wTf];





function tw=antitrigauss_full(n,alpha,beta)

omega=(beta-alpha)/2;

ab=r_trigauss(n+1,omega);
ab(end,2)=2*ab(end,2);
xwAT=gauss(n+1,ab);

xAT=xwAT(:,1); wAT=xwAT(:,2);
tAT=2*asin(sin(omega/2)*xwAT(:,1))+(beta+alpha)/2;

xwT=gauss(n,ab(1:end-1,:));
xT=xwT(:,1); wT=xwT(:,2);
tT=2*asin(sin(omega/2)*xT)+(beta+alpha)/2;
tTe=[tT; 0];
wTe=[wT; 0];

tw=[tAT wAT tTe wTe];




function tw=gauss_legendre(n,alpha,beta)

% [xGL, wGL] = legpts(N); wGL=wGL';
ab=r_jacobi(n,0,0);
xw=gauss(n,ab); xGL=xw(:,1); wGL=xw(:,2);
t=(beta+alpha)/2+(beta-alpha)*xGL/2;
w=(beta-alpha)*wGL/2;
tw=[t w];




function tw=trigauss_classic(n,alpha,beta)

% by Gaspare Da Fies and Marco Vianello, University of Padova
% 8 Nov 2011

% computes the n+1 angles and weights of a trigonometric gaussian
% quadrature formula on [alpha,beta], 0<beta-alpha<=pi

% uses the routines chebyshev.m, gauss.m from
% www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
% we suggest to put the following statements
% ab = zeros(N,2); sig = zeros(N+1,2*N);
% at the beginning of the body of chebyshev.m to speed-up execution

% input:
% n: trigonometric degree of exactness
% [alpha,beta]: angular interval, 0<beta-alpha<=pi

% output:
% tw: (n+1) x 2 array of (angles,weights)

% the formula integrates the canonical trigonometric basis with accuracy
% from about 10^(-15) (small omega) to about 10^(-13) (omega-->pi)
% up to n=300


% half-length of the angular interval
omega=(beta-alpha)/2;

if omega == pi
    tw=circle_trigquad(n,0,2*pi);
else
    
    ab=r_trigauss(n,omega);
    
    % Gaussian formula for the weight function above
    xw=gauss(n+1,ab);
    
    % angles and weights for the trigonometric gaussian formula
    tw(:,1)=2*asin(sin(omega/2)*xw(:,1))+(beta+alpha)/2;
    tw(:,2)=xw(:,2);
    
end



function ab=r_trigauss(n,omega)

% modified Chebyshev moments by recurrence
z(1)=2*omega;

z(n+1)=quadgk(@(t)cos(2*n*acos(sin(t/2)/sin(omega/2))),...
    -omega,omega,'MaxIntervalCount',50000);
temp=(2:2:2*n-1);
dl=1/4-1./(4*(temp-1));
dc=1/2-1/sin(omega/2)^2-1./(2*(temp.^2-1));
du=1/4+1./(4*(temp+1));
d=4*cos(omega/2)/sin(omega/2)./(temp.^2-1)';
d(n-1)=d(n-1)-du(n-1)*z(n+1);
z(2:n)=tridisolve(dl(2:n-1),dc(1:n-1),du(1:n-2),d(1:n-1));
mom=zeros(1,2*n+2);
mom(1:2:2*n+1)=z(1:n+1);

% normalization of the moments (monic polynomials)
k=(3:length(mom));
mom(3:end)=exp((2-k)*log(2)).*mom(3:end);

% recurrence coeffs of the monic Chebyshev polynomials
abm(:,1)=zeros(2*n+1,1);
abm(:,2)=0.25*ones(2*n+1,1); abm(1,2)=pi; abm(2,2)=0.5;

% recurrence coeffs for the monic OPS w.r.t. the weight function
% w(x)=2*sin(omega/2)/sqrt(1-sin^2(omega/2)*x^2)
% by the modified Chebyshev algorithm
[ab,normsq]=chebyshev(n+1,mom,abm);






function x = tridisolve(a,b,c,d)
%   TRIDISOLVE  Solve tridiagonal system of equations.
% From Cleve Moler's Matlab suite
% http://www.mathworks.it/moler/ncmfilelist.html

%     x = TRIDISOLVE(a,b,c,d) solves the system of linear equations
%     b(1)*x(1) + c(1)*x(2) = d(1),
%     a(j-1)*x(j-1) + b(j)*x(j) + c(j)*x(j+1) = d(j), j = 2:n-1,
%     a(n-1)*x(n-1) + b(n)*x(n) = d(n).
%
%   The algorithm does not use pivoting, so the results might
%   be inaccurate if abs(b) is much smaller than abs(a)+abs(c).
%   More robust, but slower, alternatives with pivoting are:
%     x = T\d where T = diag(a,-1) + diag(b,0) + diag(c,1)
%     x = S\d where S = spdiags([[a; 0] b [0; c]],[-1 0 1],n,n)

x = d;
n = length(x);
for j = 1:n-1
    mu = a(j)/b(j);
    b(j+1) = b(j+1) - mu*c(j);
    x(j+1) = x(j+1) - mu*x(j);
end
x(n) = x(n)/b(n);
for j = n-1:-1:1
    x(j) = (x(j)-c(j)*x(j+1))/b(j);
end






function tw=circle_trigquad(n,alpha,beta)

if nargin == 1
    alpha=0;
    beta=2*pi;
end
N=n+1;
w=(2*pi/N)*ones(N,1);
t=linspace(pi/N,2*pi-pi/N,N); t=alpha+t';

tw=[t w];





function tw=anti_trigauss(n,alpha,beta)

% by Gaspare Da Fies and Marco Vianello, University of Padova
% 8 Nov 2011

% computes the n+1 angles and weights of a trigonometric gaussian
% quadrature formula on [alpha,beta], 0<beta-alpha<=pi

% uses the routines chebyshev.m, gauss.m from
% www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
% we suggest to put the following statements
% ab = zeros(N,2); sig = zeros(N+1,2*N);
% at the beginning of the body of chebyshev.m to speed-up execution

% input:
% n: trigonometric degree of exactness
% [alpha,beta]: angular interval, 0<beta-alpha<=pi

% output:
% tw: (n+1) x 2 array of (angles,weights)

% the formula integrates the canonical trigonometric basis with accuracy
% from about 10^(-15) (small omega) to about 10^(-13) (omega-->pi)
% up to n=300


% half-length of the angular interval
omega=(beta-alpha)/2;

% We're calculating Antigaussian points, so we want to shift from
% G(n) to AG(n+1)
n = n+1;

ab=r_trigauss(n,omega);

% Only other real modification from original gaussian nodes
ab(end,2) = ab(end,2)*2;

% Gaussian formula for the weight function above
xw=gauss(n+1,ab);

% angles and weights for the trigonometric gaussian formula
tw(:,1)=2*asin(sin(omega/2)*xw(:,1))+(beta+alpha)/2;
tw(:,2)=xw(:,2);









% CHEBYSHEV Modified Chebyshev algorithm.
%
%    Given a weight function w encoded by its first 2n modified
%    moments, stored in the (row) vector mom, relative to monic 
%    polynomials defined by the (2n-1)x2 array abm of their
%    recurrence coefficients, [ab,normsq]=CHEBYSHEV(n,mom,abm)
%    generates the array ab of the first n recurrence coefficients
%    of the orthogonal polynomials for the weight function w, and 
%    the vector normsq of their squared norms. The n alpha-
%    coefficients are stored in the first column, the n beta-
%    coefficients in the second column, of the nx2 array ab. The
%    call [ab,normsq]=CHEBYSHEV(n,mom) does the same, but using the 
%    classical Chebyshev algorithm. If n is larger than the sizes
%    of mom and abm warrant, then n is reduced accordingly.
%
function [ab,normsq]=chebyshev(N,mom,abm)
if N<=0, error('N out of range'), end
if N>size(mom,2)/2, N=size(mom,2)/2; end
if nargin<3, abm=zeros(2*N-1,2); end
if N>(size(abm,1)+1)/2, N=(size(abm,1)+1)/2; end
ab(1,1)=abm(1,1)+mom(2)/mom(1); ab(1,2)=mom(1);
if N==1, normsq(1)=mom(1); return, end
sig(1,1:2*N)=0; sig(2,:)=mom(1:2*N);
for n=3:N+1
  for m=n-1:2*N-n+2
    sig(n,m)=sig(n-1,m+1)-(ab(n-2,1)-abm(m,1))*sig(n-1,m) ...
      -ab(n-2,2)*sig(n-2,m)+abm(m,2)*sig(n-1,m-1);
  end
  ab(n-1,1)=abm(n-1,1)+sig(n,n)/sig(n,n-1)-sig(n-1,n-1)/ ...
    sig(n-1,n-2);
  ab(n-1,2)=sig(n,n-1)/sig(n-1,n-2);
end
for n=1:N, normsq(n)=sig(n+1,n); end; normsq=normsq'; 













function xw=gauss(N,ab)
N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
  J(n,n-1)=sqrt(ab(n,2));
  J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];