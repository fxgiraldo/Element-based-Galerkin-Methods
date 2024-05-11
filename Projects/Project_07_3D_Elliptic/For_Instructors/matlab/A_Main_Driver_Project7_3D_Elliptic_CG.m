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
nel=1;
nop=5;    %Interpolation Order
noq=nop + 1; %Integration Order
c=1; %exact solution variable
plot_grid=1; %=0 Don't plot, =1 Plot Grid
plot_solution=0;
plot_matrices=0;
lwarp_grid=0; %0=don't, 1=do
%-------------------------Only Change These Lines------------------%

ngl=nop + 1;
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
Mmatrix = create_Mmatrix(intma,jac,wnq,psi,npoin,nelem,ngl,nq);
Lmatrix = create_Lmatrix(intma,jac,wnq,ksi_x,ksi_y,ksi_z,...
          eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,psi,dpsi,...
          npoin,nelem,ngl,nq);
%------------------------Ask Students to add these functions---------%
%------------------------Ask Students to add these functions---------%

%Compute Exact Solution
[qe,fe]=exact_solution(coord,npoin,c);

%Impose Homogeneous Dirichlet Boundary Conditions
Rvector=Mmatrix*fe;
for i=1:nboun
    ip=iboun(i);
    Lmatrix(ip,:)=0;
    Lmatrix(ip,ip)=1;
    Rvector(ip)=qe(ip);
end %i

%Solve System 
q0=Lmatrix\Rvector; 

%Compute Norm
l2_norm=norm(q0-qe,2)/norm(qe,2);
t1=cputime;
dt=t1-t0;
disp([' nop  = ' num2str(nop),' noq  = ' num2str(noq),' nel = ' num2str(nel),' npoin = ' num2str(npoin),' l2 norm = ' num2str(l2_norm), ' cpu = ' num2str(dt) ]);

%Plot Exact Solution
if (plot_solution == 1)
    h=figure;
    figure(h)
    subplot(1,2,1);
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
    
    [xi,yi]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
    qi=griddata(xe,ye,ze,qe,xi,yi,zi);
    surf(xi,yi,zi,qi);
    %colorbar('SouthOutside');
    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    ylabel('Z','FontSize',18);
    axis image
    title_text=['Exact Solution For: Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);


%     [xi,yi,zi]=meshgrid(xmin:dx:xmax,ymin:dy:ymax,zmin:dz:zmax);
%     qi=griddata(xe,ye,ze,qe,xi,yi,zi);
%     surf(xi,yi,zi,qi);
%     %colorbar('SouthOutside');
%     xlabel('X','FontSize',18);
%     ylabel('Y','FontSize',18);
%     ylabel('Z','FontSize',18);
%     axis image
%     title_text=['Exact Solution For: Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq)];
%     title([title_text],'FontSize',18);
%     set(gca, 'FontSize', 18);

%     %Plot Solution
%     subplot(1,2,2);
%     qi=griddata(xe,ye,ze,q0,xi,yi,zi);
%     surf(xi,yi,zi,qi);
%     %colorbar('SouthOutside');
%     xlabel('X','FontSize',18);
%     ylabel('Y','FontSize',18);
%     ylabel('Z','FontSize',18);
%     axis image
%     title_text=['Numerical Solution: L2 Norm = ' num2str(l2_norm)];
%     title([title_text],'FontSize',18);
%     set(gca, 'FontSize', 18);
%     %file_ps=['se_n' num2str(nelem) 'p' num2str(nop)];
%     %eval(['print ' file_ps ' -depsc']);
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
