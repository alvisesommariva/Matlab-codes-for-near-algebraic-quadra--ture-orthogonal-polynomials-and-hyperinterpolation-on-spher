
function demo_hyperinterp_01

%--------------------------------------------------------------------------
% Object:
% Demo on algebraic polynomial hyperinterpolation on spherical triangles.
% The determination of the X pointset varies with the choice of the
% parameter "pts_type".
%
% This code has been used for numerical experiments on:
%
% A. Sommariva and M. Vianello
% Numerical hyperinterpolation overspherical triangles
%
% For degrees higher than 10 it may be a little time consuming. See "cpus"
% in the table below (Example section), for additional details.
%--------------------------------------------------------------------------
% Required routines:
% 1. define_domain (attached)
% 2. define_function (attached)
% 3. cub_sphtri (external,requires other subroutines)
% 4. dPOLYFITsph (external, requires "dORTHVANDsph")
% 5. dPOLYVALsph (external, requires "dCHEBVAND")
% 6. dLEBsph (external, requires "dCHEBVAND", "dORTHVANDsph")
% 7. plot_s2 (external)
%--------------------------------------------------------------------------
% Example:
% >> demo_hyperinterp_01
% 
%  	  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 
%  	 Plot cap (sph-tri) 
%  
% 
%  	 .........................................................................
%  	 ** Domain:  0
%  	 ** Hyperinterpolation on compressed pointset 
%  	
%  	 .........................................................................
%  	 |  n  |  hypset  |  initset  |  dimpoly |  cpus   |    leb   |
%  	 .........................................................................
%  	 |   1 |       9  |    4959   |     4    | 1.1e-02 | 3.60e+00 |
%  	 |   2 |      25  |    5310   |     9    | 1.2e-02 | 6.25e+00 |
%  	 |   3 |      49  |    5673   |    16    | 2.6e-02 | 9.73e+00 |
%  	 |   4 |      81  |    6048   |    25    | 3.7e-02 | 1.30e+01 |
%  	 |   5 |     121  |    6435   |    36    | 8.4e-02 | 1.62e+01 |
%  	 |   6 |     169  |    6834   |    49    | 1.3e-01 | 1.97e+01 |
%  	 |   7 |     225  |    7245   |    64    | 2.8e-01 | 2.62e+01 |
%  	 |   8 |     289  |    7668   |    81    | 4.4e-01 | 2.59e+01 |
%  	 |   9 |     361  |    8103   |   100    | 6.4e-01 | 3.20e+01 |
%  	 |  10 |     441  |    8550   |   121    | 1.5e+00 | 3.44e+01 |
%  	 |  11 |     529  |    9009   |   144    | 2.0e+00 | 3.55e+01 |
%  	 |  12 |     625  |    9480   |   169    | 2.2e+00 | 3.98e+01 |
%  	 |  13 |     729  |    9963   |   196    | 5.0e+00 | 4.34e+01 |
%  	 |  14 |     841  |   10458   |   225    | 1.1e+01 | 4.90e+01 |
%  	 |  15 |     961  |   10965   |   256    | 1.9e+01 | 5.18e+01 |
%  	 .........................................................................
%  
% 
%  	 .........................................................................
%  	 For each function "f" and each degree we display the L2 relative weighted
%  	 hyperinterp. errors REL, evaluated w.r.t. to a high order cubature.
%  	 Next we compute RELinf, i.e. "maximum absolute error / norm(f,inf)".
%  	 .........................................................................
%  
%  
%  	 *** 1+x+y.^2+x.^2.*y+x.^4+y.^5+x.^2.*y.^2.*z.^2
% 
%  	 .....................................
%  	 |  deg  |   REL2   |  RELinf  |
%  	 .....................................
%  	 |    1  | 5.48e-02 | 9.12e-02 |
%  	 |    2  | 1.93e-02 | 5.85e-02 |
%  	 |    3  | 2.83e-03 | 1.61e-02 |
%  	 |    4  | 2.86e-04 | 1.33e-03 |
%  	 |    5  | 4.70e-05 | 2.25e-04 |
%  	 |    6  | 5.40e-16 | 2.26e-15 |
%  	 |    7  | 8.31e-16 | 7.34e-15 |
%  	 |    8  | 6.26e-16 | 5.08e-15 |
%  	 |    9  | 7.13e-16 | 5.50e-15 |
%  	 |   10  | 6.97e-16 | 4.02e-15 |
%  	 |   11  | 7.88e-16 | 8.89e-15 |
%  	 |   12  | 8.71e-16 | 1.10e-14 |
%  	 |   13  | 9.45e-16 | 6.28e-15 |
%  	 |   14  | 8.91e-16 | 7.48e-15 |
%  	 |   15  | 9.51e-16 | 1.37e-14 |
%  	 .....................................
%  
%  
%  	 *** cos(10*(x+y+z))
% 
%  	 .....................................
%  	 |  deg  |   REL2   |  RELinf  |
%  	 .....................................
%  	 |    1  | 9.11e-01 | 1.66e+00 |
%  	 |    2  | 6.66e-01 | 3.80e+00 |
%  	 |    3  | 1.33e-01 | 1.82e+00 |
%  	 |    4  | 1.23e-01 | 1.48e+00 |
%  	 |    5  | 1.34e-02 | 3.74e-01 |
%  	 |    6  | 1.03e-02 | 1.69e-01 |
%  	 |    7  | 8.59e-04 | 3.76e-02 |
%  	 |    8  | 5.00e-04 | 1.17e-02 |
%  	 |    9  | 3.33e-05 | 1.58e-03 |
%  	 |   10  | 1.53e-05 | 5.26e-04 |
%  	 |   11  | 8.57e-07 | 4.35e-05 |
%  	 |   12  | 3.32e-07 | 1.32e-05 |
%  	 |   13  | 1.66e-08 | 9.78e-07 |
%  	 |   14  | 5.29e-09 | 2.42e-07 |
%  	 |   15  | 2.41e-10 | 1.65e-08 |
%  	 .....................................
%  
%  
%  	 *** exp(-g(x,y,z)), g=@(x,y,z) (x-x0).^2+(y-y0).^2+(z-z0).^2, x0=0;y0=0;z0=1;
% 
%  	 .....................................
%  	 |  deg  |   REL2   |  RELinf  |
%  	 .....................................
%  	 |    1  | 1.12e-01 | 1.33e-01 |
%  	 |    2  | 1.22e-02 | 2.80e-02 |
%  	 |    3  | 1.37e-03 | 2.99e-03 |
%  	 |    4  | 1.13e-04 | 2.96e-04 |
%  	 |    5  | 8.35e-06 | 2.18e-05 |
%  	 |    6  | 5.43e-07 | 1.65e-06 |
%  	 |    7  | 2.97e-08 | 1.06e-07 |
%  	 |    8  | 1.51e-09 | 5.83e-09 |
%  	 |    9  | 6.99e-11 | 2.76e-10 |
%  	 |   10  | 2.86e-12 | 1.39e-11 |
%  	 |   11  | 1.10e-13 | 5.79e-13 |
%  	 |   12  | 3.98e-15 | 2.10e-14 |
%  	 |   13  | 7.42e-16 | 4.16e-15 |
%  	 |   14  | 9.02e-16 | 6.49e-15 |
%  	 |   15  | 6.63e-16 | 6.77e-15 |
%  	 .....................................
%  
%  
%  	 *** exp(-g(x,y,z)), g=@(x,y,z) (x-x0).^2+(y-y0).^2+(z-z0).^2  centroid=(5.774e-01,5.774e-01,5.774e-01)
% 
%  	 .....................................
%  	 |  deg  |   REL2   |  RELinf  |
%  	 .....................................
%  	 |    1  | 2.36e-02 | 1.24e-01 |
%  	 |    2  | 1.25e-03 | 1.11e-02 |
%  	 |    3  | 5.96e-05 | 7.54e-04 |
%  	 |    4  | 2.47e-06 | 3.81e-05 |
%  	 |    5  | 8.45e-08 | 1.86e-06 |
%  	 |    6  | 2.49e-09 | 5.18e-08 |
%  	 |    7  | 6.69e-11 | 2.25e-09 |
%  	 |    8  | 1.52e-12 | 4.40e-11 |
%  	 |    9  | 3.19e-14 | 1.24e-12 |
%  	 |   10  | 9.30e-16 | 2.74e-14 |
%  	 |   11  | 6.83e-16 | 7.55e-15 |
%  	 |   12  | 8.00e-16 | 7.33e-15 |
%  	 |   13  | 8.66e-16 | 5.16e-15 |
%  	 |   14  | 8.54e-16 | 7.44e-15 |
%  	 |   15  | 8.78e-16 | 1.13e-14 |
%  	 .....................................
%  
%  
%  	 *** ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(5/2), x0=0; y0=0; z0=1;
% 
%  	 .....................................
%  	 |  deg  |   REL2   |  RELinf  |
%  	 .....................................
%  	 |    1  | 2.52e-01 | 2.74e-01 |
%  	 |    2  | 1.38e-02 | 2.40e-02 |
%  	 |    3  | 8.51e-04 | 1.40e-03 |
%  	 |    4  | 1.23e-04 | 2.75e-04 |
%  	 |    5  | 2.57e-05 | 6.98e-05 |
%  	 |    6  | 7.60e-06 | 1.85e-05 |
%  	 |    7  | 2.45e-06 | 7.48e-06 |
%  	 |    8  | 8.11e-07 | 3.12e-06 |
%  	 |    9  | 3.20e-07 | 1.27e-06 |
%  	 |   10  | 1.35e-07 | 5.59e-07 |
%  	 |   11  | 5.54e-08 | 2.91e-07 |
%  	 |   12  | 2.44e-08 | 1.27e-07 |
%  	 |   13  | 1.15e-08 | 5.04e-08 |
%  	 |   14  | 5.34e-09 | 2.92e-08 |
%  	 |   15  | 2.60e-09 | 1.68e-08 |
%  	 .....................................
%  
%  
%  	 *** ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(5/2),  centroid=(5.774e-01,5.774e-01,5.774e-01)
% 
%  	 .....................................
%  	 |  deg  |   REL2   |  RELinf  |
%  	 .....................................
%  	 |    1  | 4.58e-01 | 6.61e-01 |
%  	 |    2  | 3.92e-02 | 7.83e-02 |
%  	 |    3  | 3.03e-03 | 7.31e-03 |
%  	 |    4  | 6.48e-04 | 1.82e-03 |
%  	 |    5  | 1.99e-04 | 6.73e-04 |
%  	 |    6  | 7.49e-05 | 2.21e-04 |
%  	 |    7  | 3.34e-05 | 1.45e-04 |
%  	 |    8  | 1.63e-05 | 6.08e-05 |
%  	 |    9  | 8.62e-06 | 3.98e-05 |
%  	 |   10  | 4.84e-06 | 2.21e-05 |
%  	 |   11  | 2.88e-06 | 1.23e-05 |
%  	 |   12  | 1.79e-06 | 9.33e-06 |
%  	 |   13  | 1.17e-06 | 6.07e-06 |
%  	 |   14  | 7.63e-07 | 3.84e-06 |
%  	 |   15  | 5.22e-07 | 3.06e-06 |
%  	 .....................................

%
% .........................................................................
% Reference paper:
% A. Sommariva and M. Vianello
% Numerical hyperinterpolation overspherical triangles
%--------------------------------------------------------------------------
% Dates:
% Written on 04/01/2021: A. Sommariva and M. Vianello.
%--------------------------------------------------------------------------

clear;

%--------------------------------------------------------------------------
% domain_example: determines the spherical triangle domain to be analysed.
%      From 0 to 7, going from large domains to tiny ones (i.r. the higher
%      "example", the smaller the domain.
%--------------------------------------------------------------------------

domain_example=0;

%--------------------------------------------------------------------------
% Function to study. The variable "function_typeV" can be:
%
% case 1, f=@(x,y,z) 1+x+y.^2+x.^2.*y+x.^4+y.^5+x.^2.*y.^2.*z.^2;
% case 2, f=@(x,y,z) cos(10*(x+y+z));
% case 3  x0=0; y0=0; z0=1;
%         g=@(x,y,z) (x-x0).^2 + (y-y0).^2 + (z-z0).^2;
%         f=@(x,y,z) exp( - g(x,y,z) );
% case 4  centroid=sum(vertices,1); centroid=centroid/norm(centroid);
%         x0=centroid(1); y0=centroid(2); z0=centroid(3);
%         g=@(x,y,z) (x-x0).^2 + (y-y0).^2 + (z-z0).^2;
%         f=@(x,y,z) exp( - g(x,y,z) );
% case 5  x0=0; y0=0; z0=1;
%         f=@(x,y,z) ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(5/2);
% case 6  centroid=sum(vertices,1); centroid=centroid/norm(centroid);
%         x0=centroid(1); y0=centroid(2); z0=centroid(3);
%         f=@(x,y,z) ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(5/2);
%
% Note that the variable "function_typeV" is a vector and several
% experiments will be performed. The error of each experiment will be
% described in a final figure with all the errors will be made.
%--------------------------------------------------------------------------

function_typeV=1:6;

%--------------------------------------------------------------------------
% Degrees in numerical experiments: can be a vector.
%--------------------------------------------------------------------------

nV=1:15;

%--------------------------------------------------------------------------
% Approximation type parameter "pts_type".
%     case 1, pts_type='Hyperinterpolation full set';
%     case 2, pts_type='Hyperinterpolation compressed set';
%
% Note: the case "2" should be used for mild values of "n", say at most 15.
%--------------------------------------------------------------------------
pts_type=2;

%--------------------------------------------------------------------------
% Plot domain and nodes: do_plot=1 (yes), do_plot=0 (no).
%--------------------------------------------------------------------------
do_plot=1;







% ........................ Main code below ................................

% ....... Apply settings to define domain, pointsets, and functions .......

% Domain
vertices=define_domain(domain_example);
P1=vertices(1,:); P2=vertices(2,:); P3=vertices(3,:);

[XYZWR]=cub_sphtri(40,P1,P2,P3); XR=XYZWR(:,1:3); wr=XYZWR(:,4);

% ........ Numerical approximation, varying the degree in "nV" ............
fprintf('\n \t ');
for k=1:length(nV)
    n=nV(k);
    fprintf('%2.0f ',n);
    
    tic;
    
    % ... extract hyperinterpolation set (notice that ade=2*n) ...
    if pts_type == 2 % compressed set
        [XYZW,XYZWC,momerr]=cub_sphtri(2*n,P1,P2,P3);
        X=XYZW(:,1:3); w=XYZW(:,4);
        pts=XYZWC(:,1:3); weights=XYZWC(:,4);
    else % full set
        [XYZW]=cub_sphtri(2*n,P1,P2,P3);
        X=XYZW(:,1:3); w=XYZW(:,4);
        pts=XYZW(:,1:3); weights=XYZW(:,4);
    end
    cpus(k)=toc;
    
    % .. testing AE_L2err hyperinterpolation error for each "f" at "deg" ..
    
    for j=1:length(function_typeV)
        
        function_type=function_typeV(j);
        
        % Function to approximate
        [g,gs]=define_function(function_type,vertices);
        
        % ... evaluate function to approximate ...
        g_pts=feval(g,pts(:,1),pts(:,2),pts(:,3));
        
        % ... determine polynomial hyperinterpolant ...
        [coeff,R,jvec,dbox] = dPOLYFITsph(n,pts,weights,g_pts,[],[]);
        
        % ... evaluate hyperinterpolant at initial pointset ...
        p_X=dPOLYVALsph(n,coeff,XR,R,jvec,dbox);
        
        % ... estimating hyperinterpolant error ...
        g_X=feval(g,XR(:,1),XR(:,2),XR(:,3));
        
        AE_L2err(k,j)=sqrt(wr'*(p_X-g_X).^2);
        
        % ... estimating filtered relative error in infinity norm ...
        % tol=10^(-12); iok=find(g_X >= tol & g_X <= 1/tol);
        AE_inf(k,j)=norm(p_X-g_X,inf);
        RE_inf(k,j)=AE_inf(k,j)/norm(g_X,inf);
        
        
    end
    
    % lebesgue constant estimate
    leb(k) = dLEBsph(n,pts,weights,X,jvec,dbox);
    
    % ... parameters for statistics ...
    card_X(k)=size(X,1);  card_pts(k)=size(pts,1);
    card_polyspace(k)=length(jvec); card_sphharm(k)=(n+1)^2;
    
end



% ............................... Plots ...................................

if do_plot == 1
    
    % ....... figure 1 (domain) .......
    
    clf;
    figure(1);
    R=norm(vertices(1,:));
    title_str='';
    plot_s2('spherical-triangle',vertices,[],pts,title_str,R);
    
    % ....... figure 2 (L2 hyp. error) .......
    
    figure(2);
    all_marks = {'o','v','*','d','x','s','.','^','+','>','<','p','h'};
    
    for j=1:length(function_typeV)
        
        % ... approximating L2 norm of function "g" ...
        function_type=function_typeV(j);
        [g,gs]=define_function(function_type,vertices);
        g_X=feval(g,XR(:,1),XR(:,2),XR(:,3));
        g_L2=sqrt(wr'*(g_X).^2);
        
        % ... computing L2 hyperinterpolation relative error ...
        RE_L2err(:,j)=AE_L2err(:,j)/g_L2;
        
        % ... plotting error ...
        
        marker_type=all_marks{j};
        semilogy(nV,RE_L2err(:,j),'LineWidth',1.5,'Marker',marker_type);
        hold on;
        
    end
    L=length(function_typeV);
    if L == 6, legend('f1','f2','f3','f4','f5','f6'); end
    hold off;
    
    
    % ....... figure 3 (hyp. error inf norm) .......
    
    figure(3);
    all_marks = {'o','v','*','d','x','s','.','^','+','>','<','p','h'};
    
    for j=1:length(function_typeV)
        
        % ... plotting error ...
        
        marker_type=all_marks{j};
        semilogy(nV,RE_inf(:,j),'LineWidth',1.5,'Marker',marker_type);
        hold on;
        
    end
    L=length(function_typeV);
    if L == 6, legend('f1','f2','f3','f4','f5','f6'); end
    hold off;
    
end


% ............................ General Statistics .........................


fprintf('\n \t .........................................................................');
fprintf('\n \t ** Domain: %2.0f',domain_example)
if pts_type == 1
    fprintf('\n \t ** Hyperinterpolation on full pointset \n \t');
else
    fprintf('\n \t ** Hyperinterpolation on compressed pointset \n \t');
end

fprintf('\n \t .........................................................................');
fprintf('\n \t |  n  |  hypset  |  initset  |  dimpoly |  cpus   |    leb   |');
fprintf('\n \t .........................................................................');
for k=1:length(nV)
    fprintf('\n \t | %3.0f |   %5.0f  |   %5.0f   | %5.0f    | %1.1e | %1.2e |',...
        nV(k),card_pts(k),card_X(k),card_polyspace(k),cpus(k),leb(k));
end
fprintf('\n \t .........................................................................');
fprintf('\n \n');


% ............................ L2 Error Statistics ........................

fprintf('\n \t .........................................................................');
fprintf('\n \t For each function "f" and each degree we display the L2 relative weighted');
fprintf('\n \t hyperinterp. errors REL, evaluated w.r.t. to a high order cubature.');
fprintf('\n \t Next we compute RELinf, i.e. "maximum absolute error / norm(f,inf)".');
fprintf('\n \t .........................................................................');
for j=1:length(function_typeV)
    function_type=function_typeV(j);
    [g,gs]=define_function(function_type,vertices);
    fprintf('\n \n \n \t *** ');disp(gs);
    fprintf('\n \t .....................................');
    fprintf('\n \t |  deg  |   REL2   |  RELinf  |');
    fprintf('\n \t .....................................');
    for k=1:length(nV)
        n=nV(k);
        fprintf('\n \t |  %3.0f  | %1.2e | %1.2e |',n,RE_L2err(k,j),...
            RE_inf(k,j));
    end
    fprintf('\n \t .....................................');
end

fprintf('\n \n');





function [f,fs]=define_function(function_type,parms)

%--------------------------------------------------------------------------
% Object:
% This routine, defines a function "f" to approximate and a string "fs" for
% possible messages to the user.
%--------------------------------------------------------------------------
% Input:
% function_type: determines the function to study.
% The first five functions has been used in the paper mentioned below.
%
% case 1, f=@(x,y,z) 1+x+y.^2+x.^2.*y+x.^4+y.^5+x.^2.*y.^2.*z.^2;
% case 2, f=@(x,y,z) cos(10*(x+y+z));
% case 3  x0=0; y0=0; z0=1;
%         g=@(x,y,z) (x-x0).^2 + (y-y0).^2 + (z-z0).^2;
%         f=@(x,y,z) exp( - g(x,y,z) );
% case 4  centroid=sum(vertices,1); centroid=centroid/norm(centroid);
%         x0=centroid(1); y0=centroid(2); z0=centroid(3);
%         g=@(x,y,z) (x-x0).^2 + (y-y0).^2 + (z-z0).^2;
%         f=@(x,y,z) exp( - g(x,y,z) );
% case 5  x0=0; y0=0; z0=1;
%         f=@(x,y,z) ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(5/2);
% case 6  centroid=sum(vertices,1); centroid=centroid/norm(centroid);
%         x0=centroid(1); y0=centroid(2); z0=centroid(3);
%         f=@(x,y,z) ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(5/2);
%--------------------------------------------------------------------------
% Output:
% f: defines a function "f" to approximate;
% fs: string with the content of the function "f".
%--------------------------------------------------------------------------
% Reference paper:
% A. Sommariva and M. Vianello
% Numerical hyperinterpolation overspherical triangles
%--------------------------------------------------------------------------

switch function_type
    case 1 % Fornberg
        f=@(x,y,z) 1+x+y.^2+x.^2.*y+x.^4+y.^5+x.^2.*y.^2.*z.^2;
        fs='1+x+y.^2+x.^2.*y+x.^4+y.^5+x.^2.*y.^2.*z.^2';
        
    case 2
        f=@(x,y,z) cos(10*(x+y+z));
        fs='cos(10*(x+y+z))';
        
    case 3 % exp and - square north pole distance
        x0=0; y0=0; z0=1;
        g=@(x,y,z) (x-x0).^2 + (y-y0).^2 + (z-z0).^2;
        f=@(x,y,z) exp( - g(x,y,z) );
        fs='exp(-g(x,y,z)), g=@(x,y,z) (x-x0).^2+(y-y0).^2+(z-z0).^2, x0=0;y0=0;z0=1;';
        
    case 4 % exp and - square centroid distance
        if nargin < 1
            vertices=[0 0 1; 1 0 0; 0 1 0];
        else
            vertices=parms;
        end
        
        % ... centroid computation ...
        centroid=sum(vertices,1); centroid=centroid/norm(centroid);
        x0=centroid(1); y0=centroid(2); z0=centroid(3);
        
        % ... function ...
        g=@(x,y,z) (x-x0).^2 + (y-y0).^2 + (z-z0).^2;
        f=@(x,y,z) exp( - g(x,y,z) );
        fs='exp(-g(x,y,z)), g=@(x,y,z) (x-x0).^2+(y-y0).^2+(z-z0).^2';
        
        % ... output string ...
        x0str=num2str(x0,'%1.3e');
        y0str=num2str(y0,'%1.3e');
        z0str=num2str(z0,'%1.3e');
        fs='exp(-g(x,y,z)), g=@(x,y,z) (x-x0).^2+(y-y0).^2+(z-z0).^2';
        fs=strcat(fs,'  centroid=(',x0str,',',y0str,',',z0str,')');
        
    case 5 % north pole distance like
        x0=0; y0=0; z0=1;
        f=@(x,y,z) ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(5/2);
        fs='((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(5/2), x0=0; y0=0; z0=1;';
        
    case 6 % centroid distance like
        if nargin < 1
            vertices=[0 0 1; 1 0 0; 0 1 0];
        else
            vertices=parms;
        end
        
        % ... centroid computation ...
        centroid=sum(vertices,1); centroid=centroid/norm(centroid);
        x0=centroid(1); y0=centroid(2); z0=centroid(3);
        
        % ... function ...
        f=@(x,y,z) ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(5/2);
        
        % ... output string ...
        x0str=num2str(x0,'%1.3e');
        y0str=num2str(y0,'%1.3e');
        z0str=num2str(z0,'%1.3e');
        fs='((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(5/2), ';
        fs=strcat(fs,'  centroid=(',x0str,',',y0str,',',z0str,')');
        
end  % end: "define_function"









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

















