%---------------------------------------------------------------------%
%This code computes the 1D Advection Equation using the 
%CG and DG methods with 3rd Order RK.
%This version constructs the Global Matrices which are good for 
%comparing CG and DG.
%Written by F.X. Giraldo on July 3, 2012
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
nelem=25; %Number of Elements
nop=4;    %Interpolation Order

kstages=3; %IRK1, IRK2, IRK3
dt=25.5e-3; %time-step, fraction of one revolution
Courant_max=0.3;
time_final=10; %final time in revolutions
nplots=50; %plotting variable - Ignore
iplot_solution=1; %Switch to Plot or Not.
iplot_matrices=0;
integration_points=1; %=1 for LGL and =2 for LG
integration_type=1; %=1 is inexact and =2 is exact
space_method_type=2; %1=CG, 2=DG, 3=HDG
cgdg_method=1; %1=separate or 2=unified
ti_method=6; 
            %1=ERK1, 
            %2=ERK2, for Ne=25, N=4, dt=6.71e-3 is stable
            %3=ERK3; for Ne=25, N=4, dt=9.76e-3 is stable
            %4=ERK4; 
            %5=ERK5;
            %6=Vermeire's (10,2); for Ne=25, N=4, dt=25.5e-3 is stable

icase=1; %case number: 1 is a Gaussian, 2 is a square wave, 3 is
         %Gaussian with source and 4 is square wave with source.
xmu=0.0; %filtering strength: 1 is full strength and 0 is no filter
ifilter=100000000; %time-step frequency that the filter is applied.
diss=1;

%Store Constants
ngl=nop + 1;
ntime=time_final/dt;
npoin_dg=ngl*nelem;
npoin_cg=nop*nelem + 1;
npoin_hdg=nelem+1;

if cgdg_method == 1 %Separate
    if space_method_type == 1 
        method_text = ['CG-Separate'];
        npoin=npoin_cg;
    elseif space_method_type == 2 %DG
        method_text = ['DG'];
        npoin=npoin_dg;
    elseif space_method_type == 3 %HDG
        method_text = ['HDG'];
        npoin=npoin_dg;
    end
elseif cgdg_method == 2 %Unified
    if space_method_type == 1 
        method_text = ['CG-Unified'];
        npoin=npoin_dg;
    elseif space_method_type == 2 %DG
        method_text = ['DG'];
        npoin=npoin_dg;
    elseif space_method_type == 3 %HDG
        method_text = ['HDG'];
        npoin=npoin_dg;
    end
end

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
main_text=[method_text ': ' integration_text];

%Compute Lagrange Polynomial and derivatives
[psi,dpsi] = lagrange_basis3(ngl,nq,xgl,xnq);
   
%Compute Filter Matrix
f = filter_init(ngl,xgl,xmu);

%Create Grid
[coord,intma_cg]=create_grid_dg(ngl,nelem,xgl);
dx=coord(2,1)-coord(1,1);
u=2;
% dt=Courant_max*dx/u;
ntime=round(time_final/dt)
%dt=time_final/ntime
Courant=u*dt/dx

%Compute Exact Solution
time=0;
qe = exact_solution_dg(coord,nelem,ngl,time,icase);
fe = source_function_dg(coord,nelem,ngl,time,icase);

%Create Local/Element Mass and Differentiation Matrices
mass = create_mass_dg(coord,nelem,ngl,nq,wnq,psi);
diff_matrix = create_diff_matrix_dg(ngl,nq,wnq,psi,dpsi);

%Form INTMA_DG
intma_dg=zeros(ngl,nelem);
ip=0;
for e=1:nelem
    for i=1:ngl
        ip=ip+1;
        intma_dg(i,e)=ip;
    end
end

%Form CG->HDG
cg_to_hdg=zeros(npoin_cg,1);
ip=0;
for e=1:nelem
    ip=ip+1;
    I=intma_cg(1,e);
    cg_to_hdg(I)=ip;
end
ip=ip+1;
I=intma_cg(ngl,nelem);
cg_to_hdg(I)=ip;

%Form Periodic BC Pointers
iperiodic_hdg=zeros(npoin_hdg,1);
for i=1:npoin_hdg
    iperiodic_hdg(i)=i;
end
iperiodic_hdg(npoin_hdg)=iperiodic_hdg(1); 
iperiodic=zeros(npoin,1);
for i=1:npoin
    iperiodic(i)=i;
end

%Store INTMA
if space_method_type == 1 && cgdg_method == 1 %CG AND Separate
    intma=intma_cg;
    iperiodic(npoin)=iperiodic(1);
elseif space_method_type == 2 || cgdg_method == 2 %DG Or Unified
    intma=intma_dg;
elseif space_method_type == 3  %HDG
    intma=intma_dg;
end

%Form Global Mass and Differentiation Matrices
Mmatrix=zeros(npoin,npoin);
Dmatrix=zeros(npoin,npoin);
Cmatrix=zeros(npoin_hdg,npoin);
Fmatrix_hdg=zeros(npoin,npoin_hdg);
Cmatrix_hdg=zeros(npoin_hdg,npoin_hdg);
if (space_method_type == 3)
    [Fmatrix,Fmatrix_hdg,Cmatrix,Cmatrix_hdg] = create_Fmatrix_hdg(intma,intma_cg,npoin,npoin_hdg,nelem,ngl,u,diss,cg_to_hdg,iperiodic_hdg);
else
    Fmatrix = create_Fmatrix_dg(intma,npoin,nelem,ngl,u,diss);
end
for e=1:nelem
    for i=1:ngl
        ip=iperiodic(intma(i,e));
        for j=1:ngl
            jp=iperiodic(intma(j,e));
            Mmatrix(ip,jp)=Mmatrix(ip,jp) + mass(i,j,e);
            Dmatrix(ip,jp)=Dmatrix(ip,jp) + u*diff_matrix(i,j);
        end
    end
end
    
%Form Long Exact Solution Vector
qexact=zeros(npoin,1);
for e=1:nelem
    for i=1:ngl
        ip=intma(i,e);
        qexact(ip)=qe(i,e);
    end
end

%Initialize State Vector
q1=qexact;
q0=qexact;
qp=qexact;
iplot=round(ntime/nplots);

%------------------------------Time Integration
[alpha,beta,stages] = construct_ERK_coefficients(ti_method);

%Matrices
% Cmatrix_hdg(npoin_hdg,npoin_hdg)=1;
DFmatrix=Dmatrix-Fmatrix;
if (space_method_type == 1 && cgdg_method == 1) %CGc
    Mmatrix(npoin,npoin)=1;
    DFmatrix=Dmatrix;
end
% Amatrix_hdg=Cmatrix_hdg - dt*alpha(stages,stages)*Cmatrix*Amatrix_inv*Fmatrix_hdg;
% Amatrix_hdg(npoin_hdg,npoin_hdg)=1;
Q=zeros(npoin,stages);
% QQ=zeros(npoin_hdg,stages);
R=zeros(npoin,stages);
for itime=1:ntime
   time=time + dt;
   
%    disp(['itime =  ',num2str(itime),' time = ', num2str(time),' courant = ', num2str(Courant)]);  
   
   %-----------------------------Solve Implicit problem using SDIRK
   if (space_method_type <= 2) %CG or DG
       Q(:,1)=q0(:);
       R(:,1)=DFmatrix*Q(:,1);
       for i=2:stages
           R_sum=zeros(npoin,1);
           for j=1:i-1
               R_sum=R_sum + alpha(i,j)*R(:,j);
           end
           RR=Mmatrix*q0 + dt*R_sum;
           Amatrix=Mmatrix;
           Q(:,i)=Amatrix\RR;
           R(:,i)=DFmatrix*Q(:,i);
       end
       R_sum=zeros(npoin,1);
       for i=1:stages
           R_sum=R_sum + beta(i)*R(:,i);
       end
       qp=q0 + dt*( Mmatrix\R_sum ); 
   
   elseif (space_method_type == 3) %HDG
%        Q(:,1)=q0(:);
%        QQ(:,1)=-Cmatrix_hdg\( Cmatrix*Q(:,1) );
%        QQ(npoin_hdg,1)=QQ(1,1);
%        R(:,1)=DFmatrix*Q(:,1) - Fmatrix_hdg*QQ(:,1);
%        for i=2:stages
%            R_sum=zeros(npoin,1);
%            for j=1:i-1
%                R_sum=R_sum + alpha(i,j)*R(:,j);
%            end
%            RR=Mmatrix*q0 + dt*R_sum;
%            QQ(:,i)=-Amatrix_hdg\( Cmatrix*Amatrix_inv*RR );
%            QQ(npoin_hdg,i)=QQ(1,i);
%            Q(:,i)=Amatrix\( RR - dt*alpha(i,i)*Fmatrix_hdg*QQ(:,i) );
%            R(:,i)=DFmatrix*Q(:,i) - Fmatrix_hdg*QQ(:,i);
% 
%        end
%        R_sum=zeros(npoin,1);
%        for i=1:stages
%            R_sum=R_sum + beta(i)*R(:,i);
%        end
%        qp=q0 + dt*( Mmatrix\R_sum );
   end
   
   %Filter Solution
   if (mod(itime,ifilter) == 0)
      rhs = apply_filter_dg(qp,f,nelem,ngl);
      qp=rhs;
   end

   %Update Q
   q0=qp;
end %itime

%Compute Exact Solution
qe = exact_solution_dg(coord,nelem,ngl,time,icase);

%Compute Norm
top=0;
bot=0;

for i=1:npoin
   top=top + (q0(i)-qexact(i))^2;
   bot=bot + qexact(i)^2;
end %i
Courant
l2_norm=sqrt( top/bot )
% npoin
% nelem
% ngl
% nq


%Compute a gridpoint solution
x_sol=zeros(npoin,1);
for ie=1:nelem
   for i=1:ngl
      ip=intma(i,ie);
      x_sol(ip)=coord(i,ie);
   end 
end

%Plot Solution
if (iplot_solution == 1)
    h=figure;
    figure(h);
    plot_handle=plot(x_sol,q0,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(x_sol,qexact,'b--');
    set(plot_handle,'LineWidth',2);

    xlabel('x','FontSize',18);
    ylabel('q(x,t)','FontSize',18);

    if (diss == 0)
       file_ps=[method_text num2str(nelem) 'p' num2str(nop)];
        legend(method_text,'Exact');	
    elseif (diss == 1)
       file_ps=[method_text num2str(nelem) 'p' num2str(nop)];
       legend(method_text,'Exact');	
    end

    title_text=[main_text ': Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm) ', T = ' num2str(time)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
end

if (iplot_matrices == 1)
% 
%     figure
%     spy(Rmatrix);
%     title_text=[' Rmatrix ' ];
%     title([title_text],'FontSize',18);
%     xlabel('Columns','FontSize',18);
%     ylabel('Rows','FontSize',18);
%     set(gca, 'FontSize', 18);
%     file_ps=['DG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) 'diss' num2str(diss) '_Rmatrix'];
%     %eval(['print ' file_ps ' -depsc']);
    
%     figure
%     spy(inv(Mmatrix)');
%     title_text=[' M^{-1} ' ];
%     title([title_text],'FontSize',18);
%     xlabel('Columns','FontSize',18);
%     ylabel('Rows','FontSize',18);
%     set(gca, 'FontSize', 18);
%     file_ps=['DG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) 'diss' num2str(diss) '_Mmatrix'];
%     %eval(['print ' file_ps ' -depsc']);
    
%     figure
%     spy(Fmatrix');
%     title_text=[' Fmatrix' ];
%     title([title_text],'FontSize',18);
%     xlabel('Columns','FontSize',18);
%     ylabel('Rows','FontSize',18);
%     set(gca, 'FontSize', 18);
%     file_ps=['DG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) 'diss' num2str(diss) '_Fmatrix'];
%     %eval(['print ' file_ps ' -depsc']);
%     
%     figure
%     spy(Fmatrix_hdg');
%     title_text=[' Fmatrix-HDG' ];
%     title([title_text],'FontSize',18);
%     xlabel('Columns','FontSize',18);
%     ylabel('Rows','FontSize',18);
%     set(gca, 'FontSize', 18);
%     file_ps=['DG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) 'diss' num2str(diss) '_Fmatrix'];
%     %eval(['print ' file_ps ' -depsc']);
%     
%     figure
%     spy(inv(Amatrix)');
%     title_text=[' Amatrix^{-1}' ];
%     title([title_text],'FontSize',18);
%     xlabel('Columns','FontSize',18);
%     ylabel('Rows','FontSize',18);
%     set(gca, 'FontSize', 18);
%     file_ps=['DG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) 'diss' num2str(diss) '_Fmatrix'];
%     %eval(['print ' file_ps ' -depsc']);
    
%     figure
%     spy(Cmatrix');
%     title_text=[' Cmatrix' ];
%     title([title_text],'FontSize',18);
%     xlabel('Columns','FontSize',18);
%     ylabel('Rows','FontSize',18);
%     set(gca, 'FontSize', 18);
%     file_ps=['DG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) 'diss' num2str(diss) '_Fmatrix'];
%     %eval(['print ' file_ps ' -depsc']);
%     
%     figure
%     spy(Cmatrix_hdg');
%     title_text=[' Cmatrix-HDG' ];
%     title([title_text],'FontSize',18);
%     xlabel('Columns','FontSize',18);
%     ylabel('Rows','FontSize',18);
%     set(gca, 'FontSize', 18);
%     file_ps=['DG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) 'diss' num2str(diss) '_Fmatrix'];
%     %eval(['print ' file_ps ' -depsc']);

    figure
    if (space_method_type <= 2)
        A=Amatrix;
    elseif (space_method_type == 3)
        A=-Amatrix_hdg;
    end
    [V,E]=eig(A);
    [m,n]=size(E);
    lambda=diag(E);
    xmin=min(real(lambda));xmax=max(real(lambda));
    ymin=min(imag(lambda));ymax=max(imag(lambda));
    xlam=real(lambda);
    ylam=imag(lambda);
    plot_handle=plot(xlam,ylam,'ro');
    title_text=[main_text ': E-values with max(Re) = ', num2str(max(real(lambda))) ];
    title([title_text],'FontSize',18);
    set(plot_handle,'LineWidth',2);
    xlabel('Re','FontSize',18);
    ylabel('Im','FontSize',18);
    xmin=min(xmin,-1);
    xmax=max(xmax,+1);
    set(gca, 'FontSize', 18);
end
