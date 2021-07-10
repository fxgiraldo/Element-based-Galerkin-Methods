%---------------------------------------------------------------------%
%This code computes the 1D Advection Equation using the 
%CG and DG methods with LSRK 4th Order 5-stage time-integrator.
%This version constructs the Global Matrices which are good for 
%comparing CG and DG and looking at the eigen-values of the RHS matrix.
%Written by F.X. Giraldo on April 22, 2021
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
%-------------------------Only Change These Lines------------------%
nelem=32; %Number of Elements
nop=4;    %Interpolation Order
integration_points=1; %=1 for LGL and =2 for LG
integration_type=1; %=1 is inexact and =2 is exact
space_method_type='cg'; %CG or DG
flux_type=2; %1=centered flux and 2=upwind

Courant_max=0.1; %dt controlled by Courant_max
time_final=1; %final time in revolutions
iplot_solution=1; %Switch to Plot or Not.
iplot_matrices=0;

icase=1; %case number: 1 is a Gaussian, 2 is a square wave, 3 is
         %Gaussian with source and 4 is square wave with source.
xmu=0.05; %filtering strength: 1 is full strength and 0 is no filter
ifilter=0; %time-step frequency that the filter is applied. 0=never, 1=every time-step
%-------------------------Only Change These Lines------------------%

%Store Constants
ngl=nop + 1;

npoin_cg=nop*nelem + 1;
npoin_dg=ngl*nelem;

%Compute Interpolation and Integration Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);
if (integration_points == 1)
    integration_text=['LGL'];
    if (integration_type == 1)
        noq=nop;
    elseif (integration_type == 2)
        noq=nop+1;
    end
    nq=noq + 1;
    [xnq,wnq]=legendre_gauss_lobatto(nq);
elseif (integration_points == 2)
    integration_text=['LG'];
    noq=nop;
    nq=noq + 1;
    [xnq,wnq]=legendre_gauss(nq);
end

%Compute Lagrange Polynomial and derivatives
[psi,dpsi] = lagrange_basis(ngl,nq,xgl,xnq);

%Create Grid
[coord_cg,coord_dg,intma_cg,intma_dg,periodicity_cg,periodicity_dg]=create_grid(ngl,nelem,npoin_cg,npoin_dg,xgl);

%Form Global Matrix and Periodic BC Pointers
if strcmp(space_method_type,'cg')
    npoin=npoin_cg;
    coord=coord_cg;
    intma=intma_cg;
    periodicity=periodicity_cg;
elseif strcmp(space_method_type,'dg')
    npoin=npoin_dg;
    coord=coord_dg;
    intma=intma_dg;
    periodicity=periodicity_dg;
end
main_text=[space_method_type ': ' integration_text];

%Compute Exact Solution
time=0;
[qe,u] = exact_solution(coord,npoin,time,icase);

%Compute Courant Number
dx=coord(2)-coord(1);
dt=Courant_max*dx/u;
ntime=round(time_final/dt);
dt=time_final/ntime;
Courant=u*dt/dx;
q0=qe;
disp(['Courant = ',num2str(Courant),' dt = ',num2str(dt),' ntime = ',num2str(ntime),' time_final = ',num2str(time_final)])

%------------------------Ask Students to add these functions---------%
% % Create Local/Element Mass and Differentiation Matrices
% Me = create_mass_matrix(intma,coord,nelem,ngl,nq,wnq,psi);
% De = create_diff_matrix(ngl,nq,wnq,psi,dpsi);
% %Form Global Mass and Differentiation Matrices
% [Mmatrix,Dmatrix] = Matrix_DSS(Me,De,u,intma,periodicity,ngl,nelem,npoin);
%------------------------Ask Students to add these functions---------%

%Apply BCs
if strcmp(space_method_type,'cg')
    Mmatrix(npoin,npoin)=1; 
    Fmatrix=zeros(npoin,npoin);
elseif strcmp(space_method_type,'dg')
    if (flux_type == 1)
        Fmatrix = Fmatrix_centered_flux(intma,nelem,npoin,ngl,u);
    elseif (flux_type == 2)
        Fmatrix = Fmatrix_upwind_flux(intma,nelem,npoin,ngl,u);
    end
end
Rmatrix=Dmatrix - Fmatrix;

%Left-Multiply by Inverse Mass Matrix
Dmatrix_hat=Mmatrix\Rmatrix;
 
%Initialize State Vector
q1=qe;
q0=qe;
qp=qe;

%Time Integration
[q0,time] = time_integration(q0,Dmatrix_hat,periodicity,time,ntime,dt);

%Compute Exact Solution
[qe,u] = exact_solution(coord,npoin,time,icase);

%Compute Norm
l2_norm=norm(q0-qe,2)/norm(q0,2);
disp(['L2 = ',num2str(l2_norm),' nop = ',num2str(nop),' nelem = ',num2str(nelem),' npoin = ',num2str(npoin)])

%Plot Solution
if (iplot_solution == 1)
    h=figure;
    figure(h);
    plot_handle=plot(coord,q0,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(coord,qe,'b--');
    set(plot_handle,'LineWidth',2);

    xlabel('x','FontSize',18);
    ylabel('q(x,t)','FontSize',18);

    title_text=[main_text ': Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm) ', T = ' num2str(time)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
end

if (iplot_matrices == 1)
    %Plot E-values
    figure
    [V,EE]=eig(Dmatrix_hat);
    [m,n]=size(EE);
    for i=1:m
        E(i)=EE(i,i);
        norm_E(i)=sqrt( conj(E(i))*E(i) );
    end
    max_norm_E=max(norm_E);
    E=E/max_norm_E;
    plot(real(E),imag(E),'ro','LineWidth',2);
    title_text=[main_text ': E-values of Dhat with max(Re) = ', num2str(max(real(E))) ];
    title([title_text],'FontSize',18);
    set(plot_handle,'LineWidth',2);
    xlabel('Re','FontSize',18);
    ylabel('Im','FontSize',18);
    set(gca, 'FontSize', 18);
    axis([ -1 0.2 -1 +1]);
        
    figure
    spy(Mmatrix);
    title_text=[' M ' ];
    title([title_text],'FontSize',18);
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
           
    figure
    spy(Fmatrix);
    title_text=[' Fmatrix ' ];
    title([title_text],'FontSize',18);
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
end