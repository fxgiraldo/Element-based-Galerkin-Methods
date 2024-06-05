%---------------------------------------------------------------------%
%This driver contains the solution to Project 7: the 3D Poisson Equation using the CG method.
%
%The approached follows Algorithm 12.18 in the book whereby the global matrices
%are constructed by combining the element matrices (Alg. 12.9 and 12.10)
%and the DSS (Alg. 12.11) to construct the global Mass and global Laplacian 
%matrices. However, here we extend to 3D.
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
%-------------------------Only Change These Lines------------------%
nel=2;
nop=5;    %Interpolation Order
integration_type=1; %1=inexact, 2=exact
c=1; %exact solution variable => Wave number in each direction
plot_grid=1; %=0 Don't plot, =1 Plot Grid
plot_solution=0; %3D Plots not working
plot_matrices=0;
lwarp_grid=1; %0=don't, 1=do
%-------------------------Only Change These Lines------------------%

ngl=nop + 1;
if (integration_type == 1)
    noq=nop;
elseif (integration_type == 2)
    noq=nop+1;
end
nq=noq + 1;

t0=cputime;
nelx=nel; nely=nel; nelz=nel;
nx=nelx*nop+1;
ny=nely*nop+1;
nz=nelz*nop+1;
npoin=nx*ny*nz;
nelem=nelx*nely*nelz;
nboun=2*(nx*nz) + 2*(ny-2)*nz + 2*(nx-2)*(ny-2);

%Compute LGL Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);

%Create Grid
[coord,intma,iboun]=create_grid(npoin,nelem,nboun,nelx,nely,nelz,ngl,xgl,plot_grid,lwarp_grid);

%Compute Legendre Cardinal functions and derivatives
[psi,dpsi,xnq,wnq] = lagrange_basis(ngl,nq,xgl);

%Compute Metric Terms
[ksi_x,ksi_y,ksi_z,eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,jac] = metrics(coord,intma,psi,dpsi,nelem,ngl,nq);

%------------------------Ask Students to add these functions---------%
%------------------------Ask Students to add these functions---------%
%Create Mass and Laplacian Matrices
% Mmatrix = create_Mmatrix(intma,jac,wnq,psi,npoin,nelem,ngl,nq);
% Lmatrix = create_Lmatrix(intma,jac,wnq,ksi_x,ksi_y,ksi_z,...
%           eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,psi,dpsi,...
%           npoin,nelem,ngl,nq);

%Or can construct them both at once:
[Mmatrix,Lmatrix] = create_Global_Matrices(intma,jac,wnq, ...
                    ksi_x,ksi_y,ksi_z,eta_x,eta_y,eta_z, ...
                    zeta_x,zeta_y,zeta_z,psi,dpsi,npoin,nelem,ngl,nq);
%------------------------Ask Students to add these functions---------%
%------------------------Ask Students to add these functions---------%

%Compute Exact Solution
[qe,fe]=exact_solution(coord,npoin,c);

%Impose Homogeneous Dirichlet Boundary Conditions
Rvector=Mmatrix*fe;
for i=1:nboun
    I=iboun(i);
    Lmatrix(I,:)=0;
    Lmatrix(I,I)=1;
    Rvector(I)=qe(I);
end %i

%Solve System 
q0=Lmatrix\Rvector; 

%Compute Norm
l2_norm=norm(q0-qe,2)/norm(qe,2);
t1=cputime;
dt=t1-t0;
disp([' nop  = ' num2str(nop),' noq  = ' num2str(noq),' nel = ' num2str(nel),' npoin = ' num2str(npoin),' l2 norm = ' num2str(l2_norm), ' cpu = ' num2str(dt) ]);

%Extrema
q0_max=max(q0(:));
q0_min=min(q0(:));
qe_max=max(qe(:));
qe_min=min(qe(:));
disp([' q0_{max}  = ' num2str(q0_max),' q0_{min}  = ' num2str(q0_min), ]);
disp([' qe_{max}  = ' num2str(qe_max),' qe_{min}  = ' num2str(qe_min), ]);

%Plot Exact Solution
if (plot_solution == 1)
    h=figure;
    figure(h)
    %subplot(1,2,1);
    xmin=min(coord(1,:)); xmax=max(coord(1,:));
    ymin=min(coord(2,:)); ymax=max(coord(2,:));
    zmin=min(coord(3,:)); zmax=max(coord(3,:));
    xe=coord(1,:);
    ye=coord(2,:);
    ze=coord(3,:);

    nx=10; ny=10; nz=10;
    dx=(xmax-xmin)/nx;
    dy=(ymax-ymin)/ny;
    dz=(zmax-zmin)/nz;
    
    [xi,yi,zi]=meshgrid(xmin:dx:xmax,ymin:dy:ymax,zmin:dz:zmax);
    qi=griddata(xe,ye,ze,qe,xi,yi,zi,'linear');

    %Extrema
    qi_max=max(qi(:));
    qi_min=min(qi(:));
    disp([' qi_{max}  = ' num2str(qi_max),' qi_{min}  = ' num2str(qi_min), ]);

    xslice=0;
    yslice=0;
    zslice=0;
    slice(xi,yi,zi,qi,xslice,yslice,zslice)
    colorbar('SouthOutside');
    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    ylabel('Z','FontSize',18);
    axis image
    title_text=['Exact Solution For: Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
end

%Plot E-values
if (plot_matrices == 1)
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
    file_ps=['CG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) '_Evalues'];
    eval(['print ' file_ps ' -depsc']);
    
    figure
    spy(Lmatrix);
    title_text=[' Lmatrix ' ];
    title([title_text],'FontSize',18);
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
    file_ps=['CG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) '_Lmatrix'];
    eval(['print ' file_ps ' -depsc']);
    
    figure
    spy(Mmatrix);
    title_text=[' Mmatrix ' ];
    title([title_text],'FontSize',18);
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
    file_ps=['CG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) '_Mmatrix'];
    eval(['print ' file_ps ' -depsc']);
end

toc
