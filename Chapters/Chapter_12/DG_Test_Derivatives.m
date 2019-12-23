%---------------------------------------------------------------------%
%This code solves the 2D Poisson Equation using Unified CG/DG methods
%with tensor product of 1D basis function with either
%Exact or Inexact Integration and Using an NPOIN based data-structure
%Written by F.X. Giraldo on 9/2014
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
% close all;

tic

%Input Data
nel=5; %Number of Elements
nop=1;    %Interpolation Order

plot_matrices=0;
integration_type=2; %=1 is inexact and =2 is exact
space_method='cgc'; %=cgc for CG continuous; 
                    %=cgd for CG discontinuous;
                    %=dg for DG
dg_method='SIP';    %=LDG or SIP
mu_constant=1e4;    %Penalty term for SIP
icase=1; %case number: 1 is a Gaussian in CW, 2 is Gaussian along X, 
         %3 is Gaussian along Y, 4 is Gaussian along Diagonal, 
         %5 is square wave in CW, 6 is square wave along X.

nelx=nel;
nely=1;
nelem=nelx*nely; %Number of Elements
ngl=nop + 1;
npts=ngl*ngl;
npoin_CG=(nop*nelx + 1)*(nop*nely + 1);
npoin_DG=npts*nelem;
nboun=2*nelx + 2*nely;
nside=2*nelem + nelx + nely;

%Compute LGL Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);

if (integration_type == 1)
    noq=nop;
    integration_text = ['Inexact'];
elseif (integration_type == 2)
    noq=nop+1;
    integration_text = ['Exact'];
end
nq=noq + 1;
main_text=[space_method ':'  integration_text];

%Compute Legendre Cardinal functions and derivatives
[psi,dpsi,xnq,wnq] = lagrange_basis(ngl,nq,xgl);

%Create CG-Storage Grid
[coord_CG,intma_CG,bsido_CG] = create_grid(npoin_CG,nelem,nboun,...
                                nelx,nely,ngl,xgl);
                            
%Create CGDG-Storage Grid
[coord,intma,DG_to_CG,npoin] = create_CGDG_Storage(space_method,npoin_CG,npoin_DG,coord_CG, ...
                                      intma_CG,nelem,ngl);
                       
%Compute Metric Terms
[ksi_x,ksi_y,eta_x,eta_y,jac] = metrics(coord,intma,psi,dpsi,wnq,nelem,ngl,nq);

%Compute Side/Edge Information
[iside,jeside] = create_side(intma_CG,bsido_CG,npoin_CG,nelem,nboun,nside,ngl);
[psideh,imapl,imapr] = create_side_dg(iside,intma_CG,nside,nelem,ngl);
[nx,ny,jac_side] = compute_normals(psideh,intma_CG,coord_CG,...
                   nside,ngl,nq,wnq,psi,dpsi);
       
%Compute Exact Solution
[qe,qe_x,qe_y,qe_xx,qe_yy] = exact_solution_test_derivatives(coord,npoin,icase);

%Create RHS Vector and LMatrix  
Mmatrix = create_Mmatrix(jac,intma,psi,npoin,nelem,ngl,nq);
Dmatrix_x=zeros(npoin,npoin);
Fmatrix_x=zeros(npoin,npoin);
for i=1:npoin
    q=zeros(npoin,1);
    q(i)=1;
    rhs = create_Dmatrix(intma,jac,ksi_x,ksi_y,eta_x,eta_y,psi,dpsi,...
          npoin,nelem,ngl,nq,q);
    Dmatrix_x(i,:)=rhs(:,1);
%     Dmatrix_y(i,:)=rhs(:,2);
    rhs = compute_flux_LDG(q,psideh,nx,ny,jac_side,psi,npoin,...
           nside,ngl,nq,imapl,imapr,intma,qe_x);          
    Fmatrix_x(i,:)=rhs(:,1);
%     Fmatrix_y(i,:)=rhs(:,2);
end
% Lmatrix=Dmatrix_x*Mmatrix_inv*Dmatrix_x + Dmatrix_y*Mmatrix_inv*Dmatrix_y;

Dhatmatrix_x=Fmatrix_x + Dmatrix_x;
Dhatmatrix_x=-Dhatmatrix_x;
% Dhatmatrix_y=Fmatrix_y + Dmatrix_y;

%Solve
q_x=Mmatrix\Dhatmatrix_x*qe_x;
% q_y=Mmatrix\Dhatmatrix_y*qe;
q0=q_x;
qe=qe_xx;

%Compute Norm
top=0;
bot=0;
for i=1:npoin
    top=top + (q0(i)-qe(i))^2;
    bot=bot + qe(i)^2;
end %e
l2_norm=sqrt( top/bot );

%Plot Solution
h=figure;
figure(h);
xmin=min(coord(:,1)); xmax=max(coord(:,1));
ymin=min(coord(:,2)); ymax=max(coord(:,2));
xe=coord(:,1);
ye=coord(:,2);
nxx=100; nyy=100;
dx=(xmax-xmin)/nxx;
dy=(ymax-ymin)/nyy;
[xi,yi]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
qi=griddata(xe,ye,q0,xi,yi,'cubic');
[cl,h]=contourf(xi,yi,qi);
%surf(xi,yi,qi);
colorbar('SouthOutside');
xlabel('X','FontSize',18);
ylabel('Y','FontSize',18);
axis image
title_text=[space_method ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm)];
title([title_text],'FontSize',18);
set(gca, 'FontSize', 18);
%file_ps=['se_n' num2str(nelem) 'p' num2str(nop)];
%eval(['print ' file_ps ' -depsc']);

h=figure;
figure(h);
xe=coord(:,1);
ye=coord(:,2);
plot(xe,q0,'ro-',xe,qe,'b+-','LineWidth',2);
xlabel('X','FontSize',18);
ylabel('Q','FontSize',18);
axis square
title_text=[space_method ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm)];
title([title_text],'FontSize',18);
set(gca, 'FontSize', 18);
legend('CGDG','Exact');

disp(['space_method = ',space_method]);
disp(['nop = ',num2str(nop),'  nelem = ',num2str(nelem) ]);
q_max=max(q0);
q_min=min(q0);
disp(['L2_Norm = ',num2str(l2_norm),' q_max = ',num2str(q_max),'  q_min = ',num2str(q_min) ]);
disp(['npoin = ',num2str(npoin),' npoin_CG = ',num2str(npoin_CG),'  npoin_DG = ',num2str(npoin_DG) ]);

%Plot E-values
if (plot_matrices == 1)
    figure
    E=eig(Dmatrix_x);
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
    
    figure
    spy(Mmatrix);
    title_text=[' Mmatrix ' ];
    title([title_text],'FontSize',18);
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
    
    figure
    spy(Dmatrix_x);
    title_text=[' Dmatrix_x ' ];
    title([title_text],'FontSize',18);
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
    
%     figure
%     spy(Fmatrix_x);
%     title_text=[' Fmatrix_x ' ];
%     title([title_text],'FontSize',18);
%     xlabel('Columns','FontSize',18);
%     ylabel('Rows','FontSize',18);
%     set(gca, 'FontSize', 18);
    
end

toc
