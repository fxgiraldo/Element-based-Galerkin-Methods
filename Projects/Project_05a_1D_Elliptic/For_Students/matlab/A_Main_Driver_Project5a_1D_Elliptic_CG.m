%---------------------------------------------------------------------%
%This code computes the 1D Poisson Equation using the CG method.
%
%This is the solution to Project 5a in Element-based Galerkin Methods
%The algorithms here are described in Algorithms 8.1 and 8.2, and DSSed
%using Alg. 5.4.
%in F.X. Giraldo's Introduction to Element-based Galerkin 
%Methods using Tensor-Product Bases: Analysis, Algorithms, and Applications.
%
%Written by F.X. Giraldo on 7/2021
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
%---------------------------------------------------------------------%
nelem=10; %number of elements
nop=4; %polynomial order
iplot_solution=1; %Switch to Plot or Not.
iplot_matrices=0;
icase=1; %case number: 1 is Sin wave with Non-zero BC, 2 is Cos wave with Zero BC
integration_type=2; %=1 is inexact and =2 is exact
%---------------------------------------------------------------------%

if integration_type == 1
    noq=nop;
elseif integration_type == 2
    noq=nop+1;
end
ngl=nop + 1;
nq=noq + 1;
npoin=nop*nelem + 1;
disp(['nop=',num2str(nop),' noq= ',num2str(nop),' nelem= ',num2str(nelem),' npoin= ',npoin]);

%Compute LGL Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);

%Compute Legendre Cardinal functions and derivatives
[psi,dpsi,xnq,wnq] = lagrange_basis(ngl,nq,xgl);

%Create Grid
[coord,intma]=create_grid(ngl,nelem,xgl);

%Compute Exact Solution
[qe,qe_x,fe] = exact_solution(coord,npoin,icase);

%Create Periodic BC Pointer
iperiodic=zeros(npoin,1);
for i=1:npoin
   iperiodic(i)=i;
end

%---------------------------------------------------------------------%
%--Students add these functions---------------------------------------%
%---------------------------------------------------------------------%
%Create Mass and Differentiation Matrices
% Mmatrix = create_Mmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,iperiodic);
% Lmatrix = create_Lmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,iperiodic);
%---------------------------------------------------------------------%

%Apply Dirichlet BC
RHSmatrix=Mmatrix*fe;
RHSmatrix(1)=qe(1);
RHSmatrix(npoin)=qe(npoin);

%Compute CG Solution
q0=-Lmatrix\RHSmatrix;

%Compute Norm
l2_norm=norm(q0-qe,2)/norm(qe,2);

%Plot Solution
if (iplot_solution == 1)
    xmin=coord(1);
    xmax=coord(npoin);
    h=figure;
    figure(h);
    plot_handle=plot(coord,q0,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(coord,qe,'b--');
    set(plot_handle,'LineWidth',2);
    title_text=['CG: Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm)];
    axis([xmin xmax -1 +1]);
    
    %Plot Elements
    for ie=1:nelem
        i1=intma(ie,1);
        i2=intma(ie,ngl);
        x1=coord(i1);
        x2=coord(i2);
        y1=q0(i1);
        y2=q0(i2);
        plot(x1,y1,'ko');
        plot(x2,y2,'ko');
    end

    title([title_text],'FontSize',18);
    xlabel('x','FontSize',18);
    ylabel('q(x,t)','FontSize',18);
    legend('CG','Exact');
    set(gca, 'FontSize', 18);
%     file_ps=['se_n' num2str(nelem) 'p' num2str(nop)];
%     eval(['print ' file_ps ' -depsc']);
end

%Plot E-values
if (iplot_matrices == 1)
    figure
    E=eig(Lmatrix);
    m=length(E);
    for i=1:m
        norm_E(i)=sqrt( conj(E(i))*E(i) );
    end
    max_norm_E=max(norm_E);
    E=E/max_norm_E;
    plot_handle=plot(real(E),imag(E),'ro');
    title_text=[' E-values with max(Re) = ', num2str(max(real(E))) ];
    title([title_text],'FontSize',18);
    set(plot_handle,'LineWidth',2);
    xlabel('Re','FontSize',18);
    ylabel('Im','FontSize',18);
    set(gca, 'FontSize', 18);
    axis([ -1 +1 -1 +1]);
    file_ps=['se_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) '_Evalues'];
    eval(['print ' file_ps ' -depsc']);
    
    figure
    spy(Lmatrix);
    title_text=[' Lmatrix ' ];
    title([title_text],'FontSize',18);
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
    file_ps=['se_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) '_Dmatrix'];
    eval(['print ' file_ps ' -depsc']);
    
    figure
    spy(Mmatrix);
    title_text=[' Mmatrix ' ];
    title([title_text],'FontSize',18);
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
    file_ps=['se_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) '_Mmatrix'];
    eval(['print ' file_ps ' -depsc']);
end

%file_ps=['se_n' num2str(nelem) 'p' num2str(nop)];
%eval(['print ' file_ps ' -depsc']);

toc
