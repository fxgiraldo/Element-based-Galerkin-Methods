%---------------------------------------------------------------------%
%This code computes the 2D Elliptic Equation using the CG method.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
nel=4;
nop=4;    %Interpolation Order
noq=nop + 1;
%noq=nop;
c=1; %exact solution variable
plot_grid=1; %=0 Don't plot, =1 Plot Grid
rotate_grid=1; %1=yes, 0=no
plot_solution=0;
    
ngl=nop + 1;
nq=noq + 1;

t0=cputime;
nelx=nel; nely=nel;
nx=nelx*nop+1;
ny=nely*nop+1;
npoin=nx*ny;
nelem=nelx*nely;
nboun=2*(nx) + 2*(ny-2);
eps=1e-8;

%Compute LGL Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);

%Compute Legendre Cardinal functions and derivatives
[psi,dpsi,xnq,wnq] = lagrange_basis(ngl,nq,xgl);

%Create Grid
[coord,intma,iboun,iperiodic]=create_grid(npoin,nelem,nboun,nelx,nely,ngl,xgl,plot_grid,rotate_grid);

%Compute Metric Terms
[ksi_x,ksi_y,eta_x,eta_y,jac] = metrics(coord,intma,psi,dpsi,wnq,nelem,ngl,nq);

%Create Mass and Laplacian Matrices
% Mmatrix = create_Mmatrix(intma,jac,psi,iperiodic,npoin,nelem,ngl,nq);
% Lmatrix = create_Lmatrix(intma,jac,ksi_x,ksi_y,eta_x,eta_y,psi,dpsi,...
%           iperiodic,npoin,nelem,ngl,nq);

%Compute Exact Solution
[qe,fe]=exact_solution(coord,npoin,c);

%Impose Homogeneous Dirichlet Boundary Conditions
% Rvector=Mmatrix*fe;
% for i=1:nboun
%     ip=iboun(i);
%     Lmatrix(ip,:)=...;
%     Lmatrix(ip,ip)=...;
%     Rvector(ip)=...;
% end %i

%Solve System 
% q0=Lmatrix\Rvector; 

%Compute Norm
% top=0;
% bot=0;
% 
% for i=1:npoin
%    top=top + (q0(i)-qe(i))^2;
%    bot=bot + qe(i)^2;
% end
% l2_norm=sqrt( top/(bot + eps) );
t1=cputime;
dt=t1-t0;
disp([' nop  = ' num2str(nop),' nel = ' num2str(nel), ' cpu = ' num2str(dt) ]);

%Plot Exact Solution
if (plot_solution == 1)
    h=figure;
    figure(h)
    xmin=min(coord(1,:)); xmax=max(coord(1,:));
    ymin=min(coord(2,:)); ymax=max(coord(2,:));
    xe=coord(1,:);
    ye=coord(2,:);
    nx=100; ny=100;
    dx=(xmax-xmin)/nx;
    dy=(ymax-ymin)/ny;
    [xi,yi]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
    qi=griddata(xe,ye,qe,xi,yi,'cubic');
    contourf(xi,yi,qi);
    %colorbar('SouthOutside');
    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    axis image
    title_text=['Exact Solution For: Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);

    %Plot Solution
    h=figure;
    figure(h);
    qi=griddata(xe,ye,q0,xi,yi,'cubic');
    contourf(xi,yi,qi);
    %colorbar('SouthOutside');
    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    axis image
    title_text=['Numerical Solution For: Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
    %file_ps=['se_n' num2str(nelem) 'p' num2str(nop)];
    %eval(['print ' file_ps ' -depsc']);
end


toc
