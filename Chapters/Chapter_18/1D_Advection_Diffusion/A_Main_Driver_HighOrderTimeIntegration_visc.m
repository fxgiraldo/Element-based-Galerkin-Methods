%---------------------------------------------------------------------%
%This code computes the 1D Advection Equation using the DG method with LGL
%with 2nd or 3rd Order RK.
%Written by F.X. Giraldo on 10/2003
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

kstages=8; %RK1, RK2, RK3, or RK4, or any arbitrary order
dt=1e-3; %time-step, fraction of one revolution
%Courant_max=0.4;
time_final=0.25; %final time in revolutions
nplots=50; %plotting variable - Number of Frames ploted
iplot_movie=0;
iplot_movie_contour=0;
iplot_solution=1; %Switch to Plot or Not.
iplot_matrices=0;
store_movie=0;
store_plot=0;
integration_points=1; %=1 for LGL and =2 for LG
integration_type=1; %=1 is inexact and =2 is exact
space_method_type='cg'; %CG or DG
% nu=5e-3; %Non-dimensional Viscosity coefficient = visc*dt/dx^(2*nlaplacian)
nu=1e-2;
nlaplacian=0;

icase=1; %case number: 1 is a Gaussian, 2 is a square wave, 3 is
         %Gaussian with source and 4 is square wave with source;
         %5 is a step wave; 6 is a Gaussian AND a Square wave; and 
         %7 is a really sharp dispersive Gaussian;
         %8 is a 2pi sine wave
xmu=0.05; %filtering strength: 1 is full strength and 0 is no filter
ifilter=0; %time-step frequency that the filter is applied.
diss=1; %=1 dissipation (Rusanov Flux) and =0 no dissipation (Central Flux)

%Store Constants
ngl=nop + 1;
npoin=nop*nelem + 1;
ntime=time_final/dt;
iplot=round(ntime/nplots);

%Turn off filter for 1st Order Polynomials or for DG
if nop == 1
    %xmu=0;
end
% if space_method_type == 'dg'
%     xmu=0;
% end

method_text = [space_method_type];

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
main_text=[method_text ':' integration_text];

%Compute Lagrange Polynomial and derivatives
[psi,dpsi] = lagrange_basis3(ngl,nq,xgl,xnq);
   
%Compute Filter Matrix
f = filter_init(ngl,xgl,xmu);

%Compute Filter Matrix at Quadrature Points
fq=psi'*f;

%Create Grid
[coord,intma]=create_grid_dg(ngl,nelem,xgl);
dx=coord(2,1)-coord(1,1);
u=2;
% dt=Courant_max*dx/u;
% ntime=round(time_final/dt)
% dt=time_final/ntime
Courant=u*dt/dx
%diff=nu*dt/dx^(2*nlaplacian)
visc=nu*dx^(2*nlaplacian)/dt %k=visc*dt/dx^2=0.0625 is stable for 
                                 %N=2 Ne=25 dt=0.01 visc=1e-2
pause(1);

%Compute Exact Solution
time=0;
qe = exact_solution_dg(coord,nelem,ngl,time,icase);
fe = source_function_dg(coord,nelem,ngl,time,icase);

%Create Periodic BC Pointer
for i=1:npoin
   iperiodic(i)=i;
end
iperiodic(npoin)=iperiodic(1);

%Create Mass Matrix
if space_method_type == 'cg'
    Mmatrix = create_mass_cg(intma,coord,npoin,nelem,ngl,nq,wnq,psi,iperiodic);
    Mmatrix_inv=inv(Mmatrix);
elseif space_method_type == 'dg'
    mass=zeros(ngl,ngl,nelem);
    mass_inv=zeros(ngl,ngl,nelem);
    mass_t=zeros(ngl,ngl);
    mass = create_mass_dg(coord,nelem,ngl,nq,wnq,psi);
    for e=1:nelem
        mass_t(:,:)=mass(:,:,e);
        mass_inv(:,:,e)=inv(mass_t(:,:));
    end %ie
end

%Initialize State Vector
qp=qe; q0=qe; q1=qe;
iframe=0;

%Initialize Arrays
rhs=zeros(ngl,nelem);
rhs_t=zeros(ngl,1);
rhs_rk4=zeros(ngl,nelem,4);

%Compute RK Time-Integration Coefficients
[a0,a1,beta] = compute_ti_coefficients(kstages);

%Compute High-Order SSP Time-Integration Coefficients
alpha = compute_ti_aux(kstages-1);

%Time Integration
for itime=1:ntime
   time=time + dt;
   
   %disp(['itime =  ',num2str(itime),' time = ', num2str(time),' courant = ', num2str(Courant)]);
      
   qq=qp;
   q_temp=zeros(ngl,nelem);
   for ik=1:kstages
      
      %Create RHS vector: INVISCID
      rhs = create_rhs(qq,coord,nelem,ngl,nq,wnq,psi,dpsi,u);
       
      %Apply Communicator
      if space_method_type == 'cg' %CG
          rhs = apply_dss(rhs,intma,Mmatrix_inv,npoin,nelem,ngl,iperiodic);
      elseif space_method_type == 'dg' %DG
          rhs = create_flux(rhs,qq,nelem,ngl,u,diss);
      
          %Multiply by Inverse Mass matrix
          for e=1:nelem
              rhs(:,e)=mass_inv(:,:,e)*rhs(:,e);
          end %e
      end
      
      %Create RHS Laplacian
      if visc > 0 && nlaplacian > 0 
          %rhs_visc=zeros(ngl,nelem);
          q_visc=qq;
          %Compute Hyperviscosity
          for i=1:nlaplacian
                rhs_visc = create_laplacian(q_visc,coord,nelem,ngl,nq,wnq,dpsi);
                if space_method_type == 'cg' %CG
                    rhs_visc = apply_dss(rhs_visc,intma,Mmatrix_inv,npoin,nelem,ngl,iperiodic);
                elseif space_method_type == 'dg' %CG
                    %Multiply by Inverse Mass matrix
                    for e=1:nelem
                         rhs_visc(:,e)=mass_inv(:,:,e)*rhs_visc(:,e);
                    end %e
                end
                q_visc=rhs_visc;
          end
          visc_total=-(-1)^nlaplacian*visc; %negative here because CREATE_LAPLACIAN forms a negative
          rhs(:,:)=rhs(:,:) + visc_total*rhs_visc(:,:);
      end
      
      %Solve System
      qp=qq + dt/2*rhs;
      if (ik < kstages)
        q_temp = q_temp + alpha(kstages+1,ik)*qq;
      elseif (ik == kstages)
        q_temp = q_temp + alpha(kstages+1,ik)*qp;
      end
      qq=qp;
   end %ik
   qp=q_temp;
 
   %Filter Solution
   if (ifilter > 0) 
       if (mod(itime,ifilter) == 0)
            if space_method_type == 'dg'
              rhs = apply_filter_dg(qp,f,nelem,ngl);
            elseif space_method_type == 'cg'
              rhs = apply_filter_cg(qp,fq,coord,nelem,ngl,nq,wnq,psi);
              rhs = apply_dss(rhs,intma,Mmatrix_inv,npoin,nelem,ngl,iperiodic);
            end
            qp=rhs;
       end
   end

   %Update Q
   q0=qp;
   
   %Plot Solution
   if (iplot_solution == 1)
       if (mod(itime,10000) == 0)
          
          h=figure;
          for ie=1:nelem
           plot(coord(:,ie),qp(:,ie),'r-','LineWidth',2);
           hold on
          end

          xlabel('x');
          ylabel('q(x,t)');
          set(gca, 'FontSize', 18);
          pause;
       end %if_mod
   end
   
   if (iplot_movie == 1)
      if (mod(itime,iplot) == 0)
           iframe=iframe + 1;
           qi_movie(:,:,iframe)=qp;
           time_movie(iframe)=time;
           qe = exact_solution_dg(coord,nelem,ngl,time,icase);
           error_movie(:,:,iframe)=abs(qp(:,:) - qe(:,:));
      end
   end 
   
end %itime

%Plot Movie
if (iplot_movie == 1)
    figure;
    for i=1:iframe
        
        for e=1:nelem
            for j=1:ngl
                x(j)=coord(j,e);
                y(j)=qi_movie(j,e,i);
            end
            plot_handle=plot(x,y,'r-','LineWidth',2);
            %set(plot_handle,'LineWidth',2);
            hold on;
        end
        %axis([-1 +1 -0.25 +1.25]);
        title_text=[main_text ' Diss = ' num2str(diss) ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time_movie(i))];
        title([title_text],'FontSize',18);
        hold on;
        
        %Plot Elements
        for e=1:nelem
            x1=coord(1,e);
            x2=coord(ngl,e);
            y1=qi_movie(1,e,i);
            y2=qi_movie(ngl,e,i);
            plot(x1,y1,'ko');
            plot(x2,y2,'ko');
        end

        title([title_text],'FontSize',18);
        xlabel('x','FontSize',18);
        ylabel('q(x,t)','FontSize',18);
        set(gca, 'FontSize', 18);
        M_i=getframe(gcf);
        M(i)=M_i;
        pause(0.1);
        hold off;
    end
    if store_movie == 1
        file_movie=[main_text '_e' num2str(nelem) '_p' num2str(nop) '_q' num2str(noq) '_diss' num2str(diss) '.avi'];
        movie2avi(M,file_movie,'fps',5);
    end
end

%Plot Contour of Space-Time
if (iplot_movie_contour == 1)
    np=nelem*ngl;
    q_sol=zeros(np,iframe);
    e_sol=zeros(np,iframe);
    x_sol=zeros(np,1);
    t_sol=zeros(iframe,1);
    ip=0;
    for ie=1:nelem
        for i=1:ngl
          ip=ip+1;
          for l=1:iframe
            q_sol(ip,l)=q_sol(ip,l) + qi_movie(i,ie,l);
            e_sol(ip,l)=e_sol(ip,l) + error_movie(i,ie,l);
          end
          x_sol(ip)=coord(i,ie);
        end 
    end
    t_sol(:)=time_movie(:);
    ntt=100; nxx=100;
    tmin=min(time_movie);
    tmax=max(time_movie);
    xmin=min(min(coord));
    xmax=max(max(coord));
    dtt=(tmax-tmin)/ntt;
    dxx=(xmax-xmin)/nxx;
    qi=zeros(nxx+1,ntt+1);
    [xi,ti]=meshgrid(xmin:dxx:xmax,tmin:dtt:tmax);
    qi=griddata(x_sol,t_sol,e_sol',xi,ti,'cubic');
    %v=[-0.05:0.05:1.05];
    emin=min(min(e_sol));
    emax=max(max(e_sol));
    v=[0.05:0.05:emax];
    figure
    [cl,h]=contour(xi,ti,qi,v);
    %clabel(cl)
    colorbar('SouthOutside');
    xlabel('X','FontSize',18);
    ylabel('T','FontSize',18);
    axis normal

end

%Compute Exact Solution
qe = exact_solution_dg(coord,nelem,ngl,time,icase);

%Compute Broken L2 Norm
top=0;
bot=0;
error=zeros(ngl,nelem);

for ie=1:nelem
   for i=1:ngl
	   top=top + (q0(i,ie)-qe(i,ie))^2;
       error(i,ie)=abs(q0(i,ie)-qe(i,ie));
       bot=bot + qe(i,ie)^2;
   end %i
end %ie
l2_norm=sqrt( top/bot );

%Plot Solution
if (iplot_solution == 1)
    
    %Compute a gridpoint solution
    np=nelem*ngl;
    q_sol=zeros(np,1);
    qe_sol=zeros(np,1);
    x_sol=zeros(np,1);
    error_sol=zeros(np,1);
    ip=0;
    for ie=1:nelem
    for i=1:ngl
          ip=ip+1;
          q_sol(ip)=q_sol(ip) + q0(i,ie);
          qe_sol(ip)=qe_sol(ip) + qe(i,ie);
          x_sol(ip)=coord(i,ie);
          error_sol(ip)=error_sol(ip) + error(i,ie);
    end 
    end
    
    h=figure;
    figure(h);
    plot_handle=plot(x_sol,q_sol,'r-','LineWidth',2);
    hold on
    plot_handle=plot(x_sol,qe_sol,'b--','LineWidth',2);
    %axis([-1 +1 -0.25 +1.25]);

    %Plot Discontinuous Elements
    for ie=1:nelem
        x1=coord(1,ie);
        x2=coord(ngl,ie);
        y1=q0(1,ie);
        y2=q0(ngl,ie);
        plot(x1,y1,'ko');
        plot(x2,y2,'ko');
    end

    %Plot Continuous Elements
    for ie=1:nelem
        i1=intma(1,ie);
        i2=intma(ngl,ie);
        x1=x_sol(i1);
        x2=x_sol(i2);
        y1=q_sol(i1);
        y2=q_sol(i2);
        plot(x1,y1,'ko');
        plot(x2,y2,'ko');
    end

    xlabel('x','FontSize',18);
    ylabel('q(x,t)','FontSize',18);

    file_ps=[main_text '_e' num2str(nelem) '_p' num2str(nop)];
    legend(main_text,'Exact');	

    title_text=[main_text ' ' 'Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm) ', T = ' num2str(time)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
    
    if store_plot == 1
        eval(['print ' file_ps ' -depsc']);
    end
    
    %Plot Error
    figure;
    plot_handle=plot(x_sol,error_sol,'r-');
    set(plot_handle,'LineWidth',2);

    xlabel('x','FontSize',18);
    ylabel('Error','FontSize',18);

    title_text=[main_text ' ' 'Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm) ', T = ' num2str(time)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
end


dt
dx
Courant
visc
q_max=max(max(q0))
q_min=min(min(q0))
l2_norm

