function xyw = gqcircsect(n,omega,r1,r2)

% by Gaspare Da Fies and Marco Vianello, University of Padova
% 8 Nov 2011

% computes the nodes and weights of a product gaussian 
% formula on a circular annular sector centered at the origin 
% with angles in [-omega,omega]

% uses the routines:
%
% r_jacobi.m, gauss.m from
% www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
%
% trigauss.m 
% http://www.math.unipd.it/~marcov/mysoft/trigauss.m

% input: 
% n: algebraic degree of exactness
% omega: half-length of the angular interval, 0<omega<=pi
% r1,r2: internal and external radius, 0<=r1<r2

% output:
% xyw: (ceil((n+2)/2) x (n+1)) x 3 array of (xnodes,ynodes,weights) 


% trigonometric gaussian formula on the arc
tw=trigauss(n,-omega,omega);

% algebraic gaussian formula on the radial segments 
ab=r_jacobi(ceil((n+2)/2),0,0);
xw=gauss(ceil((n+2)/2),ab);
xw(:,1)=xw(:,1)*(r2-r1)/2+(r2+r1)/2;
xw(:,2)=xw(:,2)*(r2-r1)/2;

% creating the polar grid
[r,theta]=meshgrid(xw(:,1),tw(:,1));
[w1,w2]=meshgrid(xw(:,2),tw(:,2));

% nodal cartesian coordinates and weights  
xyw(:,1)=r(:).*cos(theta(:));
xyw(:,2)=r(:).*sin(theta(:));
xyw(:,3)=r(:).*w1(:).*w2(:);









% R_JACOBI Recurrence coefficients for monic Jacobi polynomials.
%
%    ab=R_JACOBI(n,a,b) generates the first n recurrence 
%    coefficients for monic Jacobi polynomials with parameters 
%    a and b. These are orthogonal on [-1,1] relative to the
%    weight function w(t)=(1-t)^a(1+t)^b. The n alpha-coefficients
%    are stored in the first column, the n beta-coefficients in
%    the second column, of the nx2 array ab. The call ab=
%    R_JACOBI(n,a) is the same as ab=R_JACOBI(n,a,a) and
%    ab=R_JACOBI(n) the same as ab=R_JACOBI(n,0,0).
%
%    Supplied by Dirk Laurie, 6-22-1998; edited by Walter
%    Gautschi, 4-4-2002.
%
function ab=r_jacobi(N,a,b)
if nargin<2, a=0; end;  if nargin<3, b=a; end
if((N<=0)|(a<=-1)|(b<=-1)) error('parameter(s) out of range'), end
nu=(b-a)/(a+b+2);
mu=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
if N==1, ab=[nu mu]; return, end 
N=N-1; n=1:N; nab=2*n+a+b;
A=[nu (b^2-a^2)*ones(1,N)./(nab.*(nab+2))];
n=2:N; nab=nab(n);
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
ab=[A' [mu; B1; B']];









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






