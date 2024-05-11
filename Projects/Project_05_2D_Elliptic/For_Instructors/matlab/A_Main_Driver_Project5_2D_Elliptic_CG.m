%---------------------------------------------------------------------%
%This driver cntains the solution to Project 5: the 2D Poisson Equation using the CG method.
%
%The approached follows Algorithm 12.18 in the book whereby the global matrices
%are constructed by combining the element matrices (Alg. 12.9 and 12.10)
%and the DSS (Alg. 12.11) to construct the global Mass and global Laplacian 
%matrices.
%
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%Cleaned up by F.X. Giraldo on 5/9/2024 
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
%-------------------------Only Change These Lines------------------%
nel=8;
nop=8;    %Interpolation Order
integration_type=1; %=1 is inexact and =2 is exact
c=1; %Wave number of exact solution in each direction
plot_grid=1; %=0 Don't plot, =1 Plot Grid
lwarp_grid=0; %1=yes, 0=no
plot_solution=1;
plot_matrices=0;
%-------------------------Only Change These Lines------------------%

t0=cputime;
nelx=nel; nely=nel;
nx=nelx*nop+1;
ny=nely*nop+1;
npoin=nx*ny;
nelem=nelx*nely;
nboun=2*(nx) + 2*(ny-2);
%eps=1e-8;

%Compute LGL Points
ngl=nop + 1;
[xgl,wgl]=legendre_gauss_lobatto(ngl);
if (integration_type == 1)
    noq=nop;
elseif (integration_type == 2)
    noq=nop+1;
end
nq=noq+1;

disp([' N  = ' num2str(nop),' Q  = ' num2str(noq),' nel = ' num2str(nel),' nelem = ' num2str(nelem),' npoin = ' num2str(npoin),' nboun = ' num2str(nboun) ]);

%Compute Legendre Cardinal functions and derivatives
[psi,dpsi,xnq,wnq] = lagrange_basis(ngl,nq,xgl);

%Create Grid
[coord,intma,iboun,periodicity]=create_grid(npoin,nelem,nboun,nelx,nely,ngl,xgl,plot_grid,lwarp_grid);

%Compute Metric Terms
[ksi_x,ksi_y,eta_x,eta_y,jac] = metrics(coord,intma,psi,dpsi,nelem,ngl,nq);

%------------------------Ask Students to add these functions---------%
%------------------------Ask Students to add these functions---------%
%Can build Mass and Laplacian Matrices separetely:
%Create Mass and Laplacian Matrices
% Mmatrix = create_Mmatrix(intma,jac,wnq,psi,periodicity,npoin,nelem,ngl,nq);
% Lmatrix = create_Lmatrix(intma,jac,wnq,ksi_x,ksi_y,eta_x,eta_y,psi,dpsi,...
%           periodicity,npoin,nelem,ngl,nq);

%Or can construct them both at once:
[Mmatrix,Lmatrix] = create_Global_Matrices(intma,jac,wnq,ksi_x,ksi_y,eta_x,eta_y,psi,dpsi,...
          periodicity,npoin,nelem,ngl,nq);
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
disp([' l2 norm = ' num2str(l2_norm), ' cpu = ' num2str(dt) ]);

%Plot Exact Solution
if (plot_solution == 1)
    h=figure;
    figure(h)
    subplot(1,2,1);
    xmin=min(coord(1,:)); xmax=max(coord(1,:));
    ymin=min(coord(2,:)); ymax=max(coord(2,:));
    xe=coord(1,:);
    ye=coord(2,:);
    nx=100; ny=100;
    dx=(xmax-xmin)/nx;
    dy=(ymax-ymin)/ny;
    [xi,yi]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
    qi=griddata(xe,ye,qe,xi,yi,'cubic');
    surf(xi,yi,qi);
    %colorbar('SouthOutside');
    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    axis image
    title_text=['Exact Solution For: Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);

    %Plot Solution
    subplot(1,2,2);
    qi=griddata(xe,ye,q0,xi,yi,'cubic');
    surf(xi,yi,qi);
    %colorbar('SouthOutside');
    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    axis image
    title_text=['Numerical Solution: L2 Norm = ' num2str(l2_norm)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
    %file_ps=['se_n' num2str(nelem) 'p' num2str(nop)];
    %eval(['print ' file_ps ' -depsc']);
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
