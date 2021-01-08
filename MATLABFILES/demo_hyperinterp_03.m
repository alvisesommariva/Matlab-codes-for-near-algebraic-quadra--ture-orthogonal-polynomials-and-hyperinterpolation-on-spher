

function demo_hyperinterp_03

%--------------------------------------------------------------------------
% Object:
% Estimate of Hyperinterpolation Lebesgue constant and estimate based on
% Christoffel function and area, over spherical rectangles.
%--------------------------------------------------------------------------
% Example:
%
% * domain_example=0;
% * nV=1:15;
%
% >> demo_hyperinterp_03
% 
%  	  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 
%  	 ...................................................
%  	 | deg |   CnC    |    CnF   |    leb   |  leb_ub  |  
%  	 ...................................................
%  	 |  1  | 1.09e+01 | 1.09e+01 | 3.61e+00 | 4.15e+00 |
%  	 |  2  | 6.36e+01 | 6.36e+01 | 6.31e+00 | 9.99e+00 |
%  	 |  3  | 2.11e+02 | 2.11e+02 | 9.90e+00 | 1.82e+01 |
%  	 |  4  | 5.26e+02 | 5.26e+02 | 1.34e+01 | 2.88e+01 |
%  	 |  5  | 1.10e+03 | 1.10e+03 | 1.68e+01 | 4.15e+01 |
%  	 |  6  | 2.02e+03 | 2.02e+03 | 2.05e+01 | 5.64e+01 |
%  	 |  7  | 3.41e+03 | 3.41e+03 | 2.75e+01 | 7.31e+01 |
%  	 |  8  | 5.36e+03 | 5.36e+03 | 2.74e+01 | 9.17e+01 |
%  	 |  9  | 7.98e+03 | 7.98e+03 | 3.39e+01 | 1.12e+02 |
%  	 | 10  | 1.14e+04 | 1.14e+04 | 3.65e+01 | 1.34e+02 |
%  	 | 11  | 1.55e+04 | 1.56e+04 | 3.78e+01 | 1.57e+02 |
%  	 | 12  | 2.08e+04 | 2.08e+04 | 4.23e+01 | 1.81e+02 |
%  	 | 13  | 2.70e+04 | 2.69e+04 | 4.60e+01 | 2.05e+02 |
%  	 | 14  | 3.48e+04 | 3.39e+04 | 5.16e+01 | 2.31e+02 |
%  	 | 15  | 4.29e+04 | 4.20e+04 | 5.40e+01 | 2.57e+02 |
%  	 ...................................................
%  	 Legend:
%  	 CnC: factor based on Christ. function (comp. nodes)
%  	 CnF: factor based on Christ. function (full nodes)
%  	 leb   : Lebesgue constant (estimate)
%  	 leb_ub: Lebesgue constant upper bound based on CnC
%  	 ...................................................
%  
%  	>> 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Define domain.
%--------------------------------------------------------------------------

domain_example=0;

%--------------------------------------------------------------------------
% Degrees in numerical experiments: can be a vector.
%--------------------------------------------------------------------------

nV=1:15;

% ........................ Main code below ................................

% ....... Apply settings to define domain, pointsets, and functions .......

% Domain
vertices=define_domain(domain_example);
P1=vertices(1,:); P2=vertices(2,:); P3=vertices(3,:);

[XYZWR]=cub_sphtri(35,P1,P2,P3); XR=XYZWR(:,1:3); wr=XYZWR(:,4);

% ........ Numerical approximation, varying the degree in "nV" ............
fprintf('\n \t ');
for k=1:length(nV)
    n=nV(k);
    
    fprintf('%2.0f ',n);
    
    tic;
    
    % ... extract hyperinterpolation set (notice that ade=2*n) ...
    
    [XYZW,XYZWC,momerr]=cub_sphtri(2*n,P1,P2,P3);
    X=XYZW(:,1:3); w=XYZW(:,4);
    pts=XYZWC(:,1:3); weights=XYZWC(:,4);
    
    % ... computing orthogonal basis (Tch. set based) ...
    [U,jvec,Q,R,dbox] = dORTHVANDsph(n,pts,weights);
    
    % ... evaluating orthogonal basis at large data set ...
    [V,dbox] = dCHEBVAND(n,XR,dbox);
    QXR=(V(:,jvec)/R);
    
    CaT(k)=norm(QXR.*QXR,inf);
    
    % ... computing orthogonal basis (Large set based) ...
    [U,jvec,Q,R,dbox] = dORTHVANDsph(n,X,w);
    
    % ... evaluating orthogonal basis at large data set ...
    [V,dbox] = dCHEBVAND(n,XR,dbox);
    QXR=(V(:,jvec)/R);
    
    CaL(k)=norm(QXR.*QXR,inf);
    
    
    % ... Lebesgue constant upper bound ...
    leb_ub(k)=sqrt(CaL(k)*sum(weights));
    
    % ... Lebesgue constant ...
    leb(k) = dLEBsph(n,pts,weights,XR,jvec,dbox);
    
    % Cm(k)=norm(QXR.*QXR,1)
end




% .............................. statistics ...............................

fprintf('\n \t ...................................................');
fprintf('\n \t | deg |   CnC    |    CnF   |    leb   |  leb_ub  |  ');
   fprintf('\n \t ...................................................');
for k=1:length(nV)
    n=nV(k);
    fprintf('\n \t | %2.0f  | %1.2e | %1.2e | %1.2e | %1.2e |',...
        n,CaT(k),CaL(k),leb(k),leb_ub(k));
end
fprintf('\n \t ...................................................'); 
fprintf('\n \t Legend:');
fprintf('\n \t CnC: factor based on Christ. function (comp. nodes)');
fprintf('\n \t CnF: factor based on Christ. function (full nodes)');
fprintf('\n \t leb   : Lebesgue constant (estimate)');
fprintf('\n \t leb_ub: Lebesgue constant upper bound based on CnC')
fprintf('\n \t ...................................................'); 

fprintf('\n \n ');



% .............................. plot .....................................

plot(nV,leb,'ro-',nV,leb_ub,'b*-');
hold on;
legend('leb','leb u.b.');
hold off;







function vertices=define_domain(example)

%--------------------------------------------------------------------------
% Object:
% This routine, defines the domain to analyse.
%--------------------------------------------------------------------------
% Input:
% example: determines the spherical triangle domain to be analysed. From 0
%   to 7, going from large domains to tiny ones (the higher "example", the
%   smaller the domain.
%--------------------------------------------------------------------------
% Output:
% vertices: it is a matrix whose rows are the vertices of the spherical
%           triangles.
%--------------------------------------------------------------------------


switch example
    
    case 0
        % orthant
        % LARGE SIZE (Asia+Europe+Australia)
        vertices=[0 0 1;
            1 0 0;
            0 1 0];
        
    case 1
        % vertices of the spherical-triangle, counterclockwise.
        a=0.9; b=0.5; % MEDIUM-LARGE SIZE (North-America)
        vertices=[0 0 1;
            0 sqrt(a) sqrt(1-a);
            sqrt(b) sqrt(a/2) sqrt(1-a/2-b)];
        
    case 2
        % vertices of the spherical-triangle, counterclockwise.
        a=0.5; b=0.2; % MEDIUM-LARGE SIZE (Australia)
        vertices=[0 0 1;
            0 sqrt(a) sqrt(1-a);
            sqrt(b) sqrt(a/2) sqrt(1-a/2-b)];
        
    case 3
        % vertices of the spherical-triangle, counterclockwise.
        a=0.1; b=0.2; % MEDIUM-SMALL SIZE (India)
        vertices=[0 0 1;
            0 sqrt(a) sqrt(1-a);
            sqrt(b) sqrt(a/2) sqrt(1-a/2-b)];
        
    case 4
        % vertices of the spherical-triangle, counterclockwise.
        a=0.01; b=0.02; % SMALL SIZE (Italy)
        vertices=[0 0 1;
            0 sqrt(a) sqrt(1-a);
            sqrt(b) sqrt(a/2) sqrt(1-a/2-b)];
        
    case 5
        % vertices of the spherical-triangle, counterclockwise.
        a=0.001; b=0.002; % VERY SMALL SIZE (Calabria+Campania)
        vertices=[0 0 1;
            0 sqrt(a) sqrt(1-a);
            sqrt(b) sqrt(a/2) sqrt(1-a/2-b)];
        
    case 6
        % vertices of the spherical-triangle, counterclockwise.
        a=0.0001; b=0.0002; % TINY SIZE (Venice Province)
        vertices=[0 0 1;
            0 sqrt(a) sqrt(1-a);
            sqrt(b) sqrt(a/2) sqrt(1-a/2-b)];
        
    case 7  % medium
        theta0=pi/4; theta1=pi/3;
        vertices=[0 0 1;
            cos(theta0) 0 sin(theta0);
            0 cos(theta1) sin(theta1)];
        
    case 8 % small
        theta0=pi/2-pi/8; theta1=pi/2-pi/10;
        vertices=[0 0 1;
            cos(theta0) 0 sin(theta0);
            0 cos(theta1) sin(theta1)];
        
    case 9
        
        vertices=[    -3.189827176852499e-01                         0     9.477605318951260e-01
            3.216349554603996e-01    -4.082482904638630e-01     8.543326569664301e-01
            -2.652237775149153e-03     4.082482904638631e-01     9.128670762866395e-01];
        
    otherwise
        % vertices of the spherical-triangle, counterclockwise.
        a=0.0000002; b=0.0000002; % VERY TINY SIZE (Murano Island)
        vertices=[0 0 1;
            0 sqrt(a) sqrt(1-a);
            sqrt(b) sqrt(a/2) sqrt(1-a/2-b)];
        
end