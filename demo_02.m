
function demo_02

%--------------------------------------------------------------------------
% Object:
% In this demo we test cubature of certain test functions, on spherical
% triangle, by rules with a choosen algebraic degree of precision.
%
% In the reference paper we have used the following preferences:
%
% * degV=5:5:30;
% * ftypeV=1:5;
% * domain_type=1;
%--------------------------------------------------------------------------
% Reference paper:
% A. Sommariva, M. Vianello
% "Near-algebraic Tchakaloff-like quadrature on spherical triangles"
%--------------------------------------------------------------------------
% Data:
% Last Update: 01/01/2021 by A. Sommariva.
%--------------------------------------------------------------------------

clear;

% ......................... degrees to analyse ............................
degV=5:5:30; % degV more than 15 may be time consuming for compression

% ......................... degrees to analyse ............................
ftypeV=1:6; % choose between 1 and 20 (goto bottom of the code)

% ......................... spherical triangle ............................

domain_type=1;

switch domain_type
    case 1 % large
        P1=[0 0 1];
        P2=[1 0 0];
        P3=[0 1 0];
    case 2  % medium
        theta0=pi/4; theta1=pi/3;
        P1=[0 0 1];
        P2=[cos(theta0) 0 sin(theta0)];
        P3=[0 cos(theta1) sin(theta1)];
    case 3 % small
        theta0=pi/2-pi/8; theta1=pi/2-pi/10;
        P1=[0 0 1];
        P2=[cos(theta0) 0 sin(theta0)];
        P3=[0 cos(theta1) sin(theta1)];
end

% ......................... filename for stats ............................

s=clock;
filename=strcat('domain_type_',num2str(domain_type),'_',num2str(s(4)), ...
    num2str(s(5)),num2str(floor(10000*s(6))));
fid = fopen(strcat(filename,'.txt'),'w');

% ........................... reference rule ..............................

nref=max(max(degV)+10,30);
[xyzwR] = cub_sphtri(nref,P1,P2,P3,0);

% ........................... numerical test ..............................

% compute rule at degree "n"
for k=1:length(degV)
    deg=degV(k);
    tic;
    [xyzw,xyzwc,momerr(k)] = cub_sphtri(deg,P1,P2,P3,0);
    cpuV(k)=toc;
    cardV(k,:)=[size(xyzw,1) size(xyzwc,1)];
    
    fprintf('\n \n \t ............ deg: %3.0f  ............ ',deg);
    fprintf('\n \t Card: Full: %6.0f Comp: %6.0f',...
        cardV(k,1),cardV(k,2));
    fprintf('\n \t Momerr: %1.3e',momerr(k));
    fprintf('\n \t CPU: Full+Comp: %1.3e',cpuV(k));
    
    for j=1:length(ftypeV)
        
        % .... function to study ....
        ftype=ftypeV(j);
        [f,fstr]=test_function(ftype);
        
        % .... cubature rule (full) ....
        nodes=xyzw(:,1:3); w=xyzw(:,4);
        fx=feval(f,nodes(:,1),nodes(:,2),nodes(:,3));
        I=fx'*w;
        
        % .... cubature rule (comp.) ....
        nodesC=xyzwc(:,1:3); wC=xyzwc(:,4);
        fxC=feval(f,nodesC(:,1),nodesC(:,2),nodesC(:,3));
        IC=fxC'*wC;
        
        % .... cubature rule (ref.) ....
        nodesR=xyzwR(:,1:3); wR=xyzwR(:,4);
        fxR=feval(f,nodesR(:,1),nodesR(:,2),nodesR(:,3));
        IR=fxR'*wR;
        
        
        % .... errors ....
        
        % absolute error
        AE(k,j)=abs(I-IR); AEC(k,j)=abs(IC-IR);
        
        % relative error
        RE(k,j)=abs(I-IR)/abs(IR); REC(k,j)=abs(IC-IR)/abs(IR);
        
        % filtered relative error
        tol=10^(-12);
        if (abs(I)>tol)|(abs(I)<1/tol)
            REF(k,j)=RE(k,j); RECF(k,j)=REC(k,j);
        else
            REF(k,j)=0; RECF(k,j)=0;
        end
        
        
        % .... statistics ....
        fprintf('\n \n \n \t *** '); disp(fstr);
        fprintf('\n \t AE: Full: %1.1e Comp: %1.1e',AE(k,j),AEC(k,j));
        fprintf('\n \t RE: Full: %1.1e Comp: %1.1e',RE(k,j),REC(k,j));
        fprintf('\n \t REF: Full: %1.1e Comp: %1.1e',REF(k,j),REC(k,j));
        
    end
end

fprintf('\n \n \n');

% .... final statistics ....

% ... filename for stats ...

s=clock;
filename=strcat('example_',num2str(domain_type),'_',num2str(s(4)), ...
    num2str(s(5)),num2str(floor(10000*s(6))));
fid = fopen(strcat(filename,'.txt'),'w');


for j=1:length(ftypeV)
    ftype=ftypeV(j);
    [f,fstr]=test_function(ftype);
    fstrf=strcat(fid,'f(x,y,z)=',fstr);
    fprintf(fid,'\n \n        * '); fprintf(fid,fstrf);
    fprintf(fid,'\n \t ...........................');
    fprintf(fid,'\n \t | deg |   ref   |  refc   |');
    fprintf(fid,'\n \t ...........................');
    for k=1:length(degV)
        fprintf(fid,'\n \t | %3.0f | %1.1e | %1.1e |',degV(k),REF(k,j),...
            REC(k,j));
    end
    fprintf(fid,'\n \t ...........................');
end

fprintf(fid,'\n \n');

s=fclose(fid);




function [f,fstr]=test_function(function_type)

switch function_type
    case 1 % Fornberg
        f=@(x,y,z) 1+x+y.^2+x.^2.*y+x.^4+y.^5+x.^2.*y.^2.*z.^2;
        fstr='1+x+y.^2+x.^2.*y+x.^4+y.^5+x.^2.*y.^2.*z.^2';
    case 2 % Fornberg
        f=@(x,y,z) 0.75*exp(-( (9*x-2).^2 + (9*y-2).^2 + (9*z-2).^2 )/4)+...
            0.75*exp( -(1/49)*((9*x+1).^2) - ((9*y+1)/10) - ((9*z+1)/10) )+...
            0.5*exp(-( (9*x-7).^2 + (9*y-3).^2 + (9*z-5).^2 )/4)+...
            -0.2*exp(-( (9*x-4).^2 + (9*y-7).^2 + (9*z-5).^2 ));
        fstr='franke type';
    case 3 % Fornberg variant
        f=@(x,y,z) (1/9)*(1+tanh(9*x-9*y+9*z));
        fstr='(1/9)*(1+tanh(9*x-9*y+9*z))';
    case 4 % Fornberg variant
        f=@(x,y,z) (1/9)*(1+sign(9*x-9*y+9*z));
        fstr='(1/9)*(1+sign(9*x-9*y+9*z))';
    case 5
        f=@(x,y,z) cos(10*(x+y+z));
        fstr='cos(10*(x+y+z))';
    case 6
        f=@(x,y,z) cos(x+y+z);
        fstr='cos(x+y+z)';
    case 7
        f=@(x,y,z) x+y+z;
        fstr='x+y+z';
    case 8
        f=@(x,y,z) exp(x);
        fstr='exp(x)';
    case 9
        f=@(x,y,z) exp(x+y+z);
        fstr='exp(x+y+z)';
    case 10
        f=@(x,y,z) -5*sin(1+10*z);
        fstr='-5*sin(1+10*z)';
    case 11
        f=@(x,y,z) 1./(101-100*z);
        fstr='1./(101-100*z)';
    case 12
        f=@(x,y,z) (abs(x)+abs(y)+abs(z))/10;
        fstr='(abs(x)+abs(y)+abs(z))/10';
    case 13
        f=@(x,y,z) 1./(abs(x)+abs(y)+abs(z)+eps);
        fstr='1./(abs(x)+abs(y)+abs(z)+eps)';
    case 14
        f=@(x,y,z) sin(1+abs(x)+abs(y)+abs(z)).^2/10;
        fstr='sin(1+abs(x)+abs(y)+abs(z)).^2/10';
    case 15
        f=@(x,y,z) ones(size(x));
        fstr='ones(size(x))';
    case 16
        f=@(x,y,z) 0.0001*abs(6-x.^2-y.^2-z.^2).^2;
        fstr='0.0001*abs(6-x.^2-y.^2-z.^2).^2';
    case 17
        f=@(x,y,z) peaks(x,y+z);
        fstr='peaks(x,y+z)';
    case 18
        f=@(x,y,z)(1.25+cos(5.4*y)).*cos(6*z)./( 6+6*((3*x-1).^2) );
        fstr='(1.25+cos(5.4*y)).*cos(6*z)./( 6+6*((3*x-1).^2) )';
    case 19
        f=@(x,y,z) exp( -(81/16)*((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2) )/3;
        fstr='exp( -(81/16)*((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2) )/3';
    case 20
        f=@(x,y,z) exp( -(81/4)*((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2))/3;
        fstr='exp( -(81/4)*((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2))/3';
    case 21
        f=@(x,y,z) x.^2+y.^2+z.^2;
        fstr='x.^2+y.^2+z.^2';
    case 22
        f=@(x,y,z) exp(-x.^2-y.^2-z.^2);
        fstr='exp(-x.^2-y.^2-z.^2)';
    case 24 % Fornberg
        f=@(x,y,z) (1/9)*(1+tanh(9*z-9*x-9*y));
        fstr='(1/9)*(1+tanh(9*z-9*x-9*y))';
    case 25 % Fornberg
        f=@(x,y,z) (1/9)*(1+sign(9*z-9*x-9*y));
        fstr='(1/9)*(1+sign(9*z-9*x-9*y))';
    case 26
        f=@(x,y,z) exp(-x.^2-5*y.^2+0.2*z.^2);
        fstr='exp(-x.^2-5*y.^2+0.2*z.^2)';
    case 27
        f=@(x,y,z) franke(x,y+z);
        fstr='franke(x,y+z)';
end