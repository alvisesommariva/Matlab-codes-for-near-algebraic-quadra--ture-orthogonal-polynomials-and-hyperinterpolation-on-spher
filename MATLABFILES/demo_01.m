
function demo_01

%--------------------------------------------------------------------------
% Object:
% In this demo we test cubature over spherical harmonics, on a spherical
% triangle, by rules with a choosen algebraic degree of precision.
%
% In the reference paper we have used the following preferences:
% * degV=5:5:30;
% * "example=1" or "example=2" or "example=3".
%
% The procedure may save summary results on files and figures.
%--------------------------------------------------------------------------
% Example:
% >> demo_01
% 
%  
%  	 ............ deg:   5  ............ 
%  	 AE Mean: Full: 3.2e-16 Comp: 3.4e-16
%  	 RE Mean filtered: Full: 3.2e-15 Comp: 3.5e-15
%  	 Card: Full:   5580 Comp:     36
%  	 Compression error: 1.500e-15
%  	 CPU: Full+Comp: 1.803e-02
%  	 Note:  11 integrals of  36 were too small or large
%  
%  	 ............ deg:  10  ............ 
%  	 AE Mean: Full: 3.2e-16 Comp: 4.4e-16
%  	 RE Mean filtered: Full: 9.5e-15 Comp: 1.4e-14
%  	 Card: Full:   6435 Comp:    121
%  	 Compression error: 8.274e-15
%  	 CPU: Full+Comp: 8.695e-02
%  	 Note:  41 integrals of 121 were too small or large
%  
%  	 ............ deg:  15  ............ 
%  	 AE Mean: Full: 1.9e-16 Comp: 3.2e-16
%  	 RE Mean filtered: Full: 1.4e-14 Comp: 2.9e-14
%  	 Card: Full:   7560 Comp:    256
%  	 Compression error: 1.771e-14
%  	 CPU: Full+Comp: 2.792e-01
%  	 Note:  88 integrals of 256 were too small or large
%  
%  	 ............ deg:  20  ............ 
%  	 AE Mean: Full: 1.3e-16 Comp: 3.9e-16
%  	 RE Mean filtered: Full: 1.1e-14 Comp: 4.6e-14
%  	 Card: Full:   8550 Comp:    441
%  	 Compression error: 2.613e-14
%  	 CPU: Full+Comp: 1.458e+00
%  	 Note: 156 integrals of 441 were too small or large
%  
%  	 ............ deg:  25  ............ 
%  	 AE Mean: Full: 1.3e-16 Comp: 5.0e-16
%  	 RE Mean filtered: Full: 1.4e-14 Comp: 5.9e-14
%  	 Card: Full:   9840 Comp:    676
%  	 Compression error: 4.716e-14
%  	 CPU: Full+Comp: 4.089e+00
%  	 Note: 241 integrals of 676 were too small or large
%  
%  	 ............ deg:  30  ............ 
%  	 AE Mean: Full: 1.7e-16 Comp: 7.6e-16
%  	 RE Mean filtered: Full: 2.8e-14 Comp: 1.2e-13
%  	 Card: Full:  10965 Comp:    961
%  	 Compression error: 8.825e-14
%  	 CPU: Full+Comp: 2.794e+01
%  	 Note: 346 integrals of 961 were too small or large
% 
%--------------------------------------------------------------------------
% Reference paper:
% A. Sommariva, M. Vianello
% "Near-algebraic Tchakaloff-like quadrature on spherical triangles"
%--------------------------------------------------------------------------
% Data:
% Last Update: 01/01/2021 by A. Sommariva.
%--------------------------------------------------------------------------


clf; clear;

% ......................... degrees to analyse ............................
degV=5:5:30;

% ......................... spherical triangle ............................

example=1;

switch example
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

% ...................  save results as files and figures...................

save_results=0; % save_results: 0 (no save), save_results: 1 (save)

% ......................... filename for stats ............................

if save_results
    s=clock;
    filename=strcat('example_',num2str(example),'_',num2str(s(4)), ...
        num2str(s(5)),num2str(floor(10000*s(6))));
    fid = fopen(strcat(filename,'.txt'),'w');
end
% ........................... reference rule ..............................

nref=max(degV)+10;
[xyzwR] = cub_sphtri(nref,P1,P2,P3,0);


% ........................... numerical test ..............................

% compute rule at degree "n"
for k=1:length(degV)
    deg=degV(k);
    tic;
    [xyzw,xyzwc,momerr(k)] = cub_sphtri(deg,P1,P2,P3,0);
    cpuV(k)=toc;
    
    % .... cubature rule (full) ....
    nodes=xyzw(:,1:3); w=xyzw(:,4);
    V=sphharmVAND(deg,nodes);
    I=V'*w;
    
    % .... cubature rule (comp.) ....
    nodesC=xyzwc(:,1:3); wC=xyzwc(:,4);
    VC=sphharmVAND(deg,nodesC);
    IC=VC'*wC;
    
    % .... cubature rule (ref.) ....
    nodesR=xyzwR(:,1:3); wR=xyzwR(:,4);
    VR=sphharmVAND(deg,xyzwR);
    IR=VR'*wR;
    
    
    % .... errors ....
    AEmean(k,:)=[mean(abs(I-IR)) mean(abs(IC-IR))];
    
    % filtered relative errors
    tol=10^(-12); iok=find(abs(I) > tol & abs(I) < 1/tol);
    If=I(iok); IRf=IR(iok); ICf=IC(iok);
    REmeanfilt(k,:)=[mean(abs(If-IRf)./abs(IRf)) ...
        mean(abs(ICf-IRf)./abs(IRf))];
    
    % .... cardinalities ....
    cardV(k,:)=[length(w) length(wC)];
    
    % .... statistics ....
    
    fprintf('\n \n \t ............ deg: %3.0f  ............ ',deg);
    fprintf('\n \t AE Mean: Full: %1.1e Comp: %1.1e',...
        AEmean(k,1),AEmean(k,2));
    fprintf('\n \t RE Mean filtered: Full: %1.1e Comp: %1.1e',...
        REmeanfilt(k,1),REmeanfilt(k,2));
    fprintf('\n \t Card: Full: %6.0f Comp: %6.0f',...
        cardV(k,1),cardV(k,2));
    fprintf('\n \t Compression error: %1.3e',momerr(k));
    fprintf('\n \t CPU: Full+Comp: %1.3e',cpuV(k));
    if length(I)-length(If) > 0
        dofilt(k)=1;
        fprintf('\n \t Note: %3.0f integrals of %3.0f were too small or large',...
            length(I)-length(If),length(I));
    else
        dofilt(k)=0;
    end
end

% .... summary ....

if save_results
    
    fprintf(fid,'\n \n \n \t');
    fprintf(fid,'                              Summary ');
    fprintf(fid,'\n \t ....................................................................');
    fprintf(fid,'\n \t |  DEG  |   REF   |   REC    |   CARD  |  CARDC  |  FILT |');
    fprintf(fid,'\n \t ....................................................................');
    for k=1:length(degV)
        
        if dofilt(k) == 1
            fprintf(fid,'\n \t |  %3.0f  | %1.1e |  %1.1e |  %6.0f | %6.0f  |   *   |',...
                degV(k),REmeanfilt(k,1),REmeanfilt(k,2),cardV(k,1),cardV(k,2));
        else
            fprintf(fid,'\n \t |  %3.0f  | %1.1e |  %1.1e |  %6.0f | %6.0f  |       |',...
                degV(k),REmeanfilt(k,1),REmeanfilt(k,2),cardV(k,1),cardV(k,2));
        end
    end
    
    fprintf(fid,'\n \t ....................................................................');
    fprintf(fid,'\n \t DEG: Degree of precision of the rule.');
    fprintf(fid,'\n \t REF: Max filtered relative error of the full rule.');
    fprintf(fid,'\n \t REC: Max filtered relative error of the compressed rule.');
    fprintf(fid,'\n \t CARD: Cardinality of the full rule.');
    fprintf(fid,'\n \t CARDC: Cardinality of the compressed rule.');
    fprintf(fid,'\n \t FILT: "*" means that some reference integrals were out of scale.');
    fprintf(fid,'\n \t ....................................................................');
    
    fclose(fid);
    
    fprintf('\n \n');
    
end
% .... plot figures ....

figure(1)
semilogy(degV,REmeanfilt(:,1),'b-','LineWidth',2); hold on;
semilogy(degV,REmeanfilt(:,2),'r-','LineWidth',2);
legend('Full rule','Comp. rule');
title('Relative filtered error');
hold off;
if save_results, saveas(gcf, strcat(filename,'A'), 'epsc'); end

figure(2)
plot(degV,cardV(:,1),'b-','LineWidth',2); hold on;
plot(degV,cardV(:,2),'r-','LineWidth',2); hold off;
legend('Full rule','Comp. rule');
title('Cardinality');
if save_results, saveas(gcf, strcat(filename,'B'), 'epsc'); end

% figure(4)
% plot_s2('spherical-triangle',[P1; P2; P3],nodes,nodesC);