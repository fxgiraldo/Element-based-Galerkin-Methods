%---------------------------------------------------------------------%
%This code computes the 1D Shallow Water Equations using either CG or 
%DG method using Local Element-Wise (LEW) Storage with either 
%Legendre or Lobatto points for interpolation and integration.
%The time-integration is accomplished via 2nd, 3rd Order, or 3rd Order 4-stage
%RK.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
nelem=8; %Number of Elements
nop=8;    %Interpolation Order

kstages=4; %RK1, RK2, RK3, or RK4
dt=1e-4; %time-step, fraction of one revolution
time_final=1.0; %final time in revolutions
nplots=25; %plotting variable - Number of Frames ploted
iplot_movie=1; %1=spit out h, U, Error, and Mass loss
               %2=spit out Error and Mass Loss separately (for book)
store_movie=0;
integration_points=1; %=1 for LGL and =2 for LG
integration_type=2; %=1 is inexact and =2 is exact
space_method_type=2; %=1 for CG and =2 for DG

icase=4; %case number: 1 is a Gaussian with flat bottom; 2 is Gaussian with linear bottom
         %3 is Gaussian with Parabolic bottom; 4 is Standing Wave (linear)
         %with Analytic Solution
xmu=0.05; %filtering strength: 1 is full strength and 0 is no filter
ifilter=0; %time-step frequency that the filter is applied.
filter_type=2; %=1 is Modal Hierarchical and =2 is regular Legendre
diss=1; %=1 dissipation (Rusanov Flux) and =0 no dissipation (Central Flux)
delta_nl=0; %=0 linear and =1 nonlinear

if (icase == 4)
    %delta_nl=0;
end

%Store Constants
eps=1e-15;
ngl=nop + 1;
npoin=nop*nelem + 1;
ntime=time_final/dt;
iplot=round(ntime/nplots);

%Turn off filter for 1st Order Polynomials
if nop == 1
    %xmu=0;
end

%Set Plotting Text for Spatial Method
if space_method_type == 1
    method_text = ['CG'];
else
    method_text = ['DG'];
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
main_text=[method_text ':' integration_text];

%Compute Lagrange Polynomial and derivatives
[psi,dpsi] = lagrange_basis3(ngl,nq,xgl,xnq);
   
%Compute Filter Matrix
f = filter_init(ngl,xgl,xmu,filter_type);

%Create Grid
[coord,intma]=create_grid_dg(ngl,nelem,xgl,icase);
dx_min=coord(2,1)-coord(1,1);

%Compute Exact Solution
time=0;
qe=zeros(2,ngl,nelem);
qp=zeros(2,ngl,nelem);
q1=zeros(2,ngl,nelem);
q0=zeros(2,ngl,nelem);
qb=zeros(ngl,nelem);
[qe,qb,gravity] = exact_solution_dg(coord,nelem,ngl,time,icase,eps);
        
%Compute Initial Mass
m1=0; m2=0;
for e=1:nelem
   dx=coord(ngl,e)-coord(1,e);
   jac=dx/2;
    for l=1:nq
      wq=wnq(l)*jac;
      for j=1:ngl
         h=qe(1,j,e)+qb(j,e);
         U=qe(2,j,e);
         m1=m1 + wq*h*psi(j,l);
         m2=m2 + wq*(U)*psi(j,l);
      end
    end
end
mass0=m1;
energy0=m2;

%Create Periodic BC Pointer
for i=1:npoin
   iperiodic(i)=i;
end
%iperiodic(npoin)=iperiodic(1);

%Create Mass Matrix
if space_method_type == 1
    [Mmatrix] = create_mass_cg(intma,coord,npoin,nelem,ngl,nq,wnq,psi,iperiodic);
    Mmatrix_inv=inv(Mmatrix);
elseif space_method_type == 2
    mass=zeros(ngl,ngl,nelem);
    mass_inv=zeros(ngl,ngl,nelem);
    mass_t=zeros(ngl,ngl);
    mass = create_mass_dg(coord,nelem,ngl,nq,wnq,psi);
    for e=1:nelem
       mass_t(:,:)=mass(:,:,e);
       mass_inv(:,:,e)=inv(mass_t(:,:));
    end %e
end

%Initialize State Vector
qp=qe; q0=qe; q1=qe;
iframe=0;

%Initialize Arrays
rhs=zeros(2,ngl,nelem);
rhs_t=zeros(ngl,1);
rhs_rk4=zeros(2,ngl,nelem,4);

%Compute RK Time-Integration Coefficients
[a0,a1,beta] = compute_ti_coefficients(kstages);

%Time Integration
for itime=1:ntime
   time=time + dt;
   u_max=-100000;
   for e=1:nelem
       for i=1:ngl
           u_max=max(u_max, qp(2,i,e)/(qp(1,i,e)+qb(i,e)));
       end
   end
   courant=u_max*dt/dx_min;
    
   if mod(itime,100) == 0
    disp(['itime =  ',num2str(itime),' time = ', num2str(time),' courant = ', num2str(courant)]);
   end
   
   for ik=1:kstages
      
      %Create RHS Matrix
      rhs = create_rhs_volume(qp,qb,coord,nelem,ngl,nq,wnq,psi,dpsi,gravity,delta_nl);
      
      %Apply Communicator
      if space_method_type == 1 %CG
          rhs = apply_bcs_cg(rhs,qp,qb,gravity,nelem,ngl,delta_nl);
          rhs = apply_dss(rhs,intma,Mmatrix_inv,npoin,nelem,ngl,iperiodic);
      elseif space_method_type == 2 %DG
          rhs = create_rhs_flux(rhs,qp,qb,nelem,ngl,diss,gravity,delta_nl);
          
          %%Multiply by Inverse Mass matrix
          for e=1:nelem
             for j=1:ngl
                for i=1:ngl
                   mass_t(i,j)=mass_inv(i,j,e);
                end %i
             end %j
             rhs_t=rhs(1,:,e);
             rhs(1,:,e)=mass_t*rhs_t(:);
             rhs_t=rhs(2,:,e);
             rhs(2,:,e)=mass_t*rhs_t(:);
          end
      end
      
      %Solve System
      qp=a0(ik)*q0 + a1(ik)*q1 + dt*beta(ik)*rhs;
                  
      %Update
      rhs_rk4(:,:,:,ik)=rhs(:,:,:);
      q1=qp;
   end %ik
   
   %Do RK4 Solution
   if (kstages == 4)
       qp(:,:,:)=q0(:,:,:) + dt/6.0*( rhs_rk4(:,:,:,1) + 2*rhs_rk4(:,:,:,2)...
                      + 2*rhs_rk4(:,:,:,3) + rhs_rk4(:,:,:,4) );
   end
   
   %Filter Solution
   if (mod(itime,ifilter) == 0)
      qp = apply_filter_dg(qp,f,nelem,ngl);
   end

   %Update Q
   q0=qp;
   
   if (iplot_movie > 0)
      if (mod(itime,iplot) == 0)
           iframe=iframe + 1;
           
           for e=1:nelem
               for i=1:ngl
                   p_movie(i,e,iframe)=qp(1,i,e);
                   u_movie(i,e,iframe)=qp(2,i,e);
               end
           end
           time_movie(iframe)=time;
           
           %Compute Mass
            m1=0;
            for e=1:nelem
               dx=coord(ngl,e)-coord(1,e);
               jac=dx/2;
                for l=1:nq
                  wq=wnq(l)*jac;
                  for j=1:ngl
                     h=qp(1,j,e)+qb(j,e);
                     U=qp(2,j,e);
                     m1=m1 + wq*h*psi(j,l);
                     m2=m2 + wq*(U)*psi(j,l);
                  end
                end
            end
            mass_movie(iframe)=abs(m1-mass0);%/mass0;
            energy_movie(iframe)=abs(m2-energy0);
            
            %Compute Norm
            [qe,qb,gravity] = exact_solution_dg(coord,nelem,ngl,time,icase,eps);
            h_top=0; h_bot=0;
            u_top=0; u_bot=0;
            for e=1:nelem
               for i=1:ngl
                   h_top=h_top + (q0(1,i,e)-qe(1,i,e))^2;
                   h_bot=h_bot + qe(1,i,e)^2;
                   u_top=u_top + (q0(2,i,e)-qe(2,i,e))^2;
                   u_bot=u_bot + qe(2,i,e)^2;
               end %i
            end %ie
            L2_norm(1,iframe)=sqrt( h_top );   
            L2_norm(2,iframe)=sqrt( u_top );   
      end
   end 
   
end %itime

if (iplot_movie == 1)
    
    figure;
    for i=1:iframe
        subplot(2,2,1); %H Solution
        hmax=max(max(max(p_movie)));
        hmin=min(min(-qb));
        xmax=max(max(coord));
        xmin=min(min(coord));
        for e=1:nelem
            for j=1:ngl
                x(j)=coord(j,e);
                y(j)=p_movie(j,e,i);
                yb(j)=-qb(j,e);
            end
            plot_handle=plot(x,y,'r-','LineWidth',2);
            hold on;
            plot_handle=plot(x,yb,'b-','LineWidth',2);
        end
        axis([xmin xmax hmin hmax]);
        title_text=[main_text ' Diss = ' num2str(diss) ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time_movie(i))];
        title([title_text],'FontSize',18);
        hold on;
         
        %Plot Elements
        for e=1:nelem
            x1=coord(1,e);
            x2=coord(ngl,e);
            y1=p_movie(1,e,i);
            y2=p_movie(ngl,e,i);
            plot(x1,y1,'ko');
            plot(x2,y2,'ko');
        end

        title([title_text],'FontSize',18);
        xlabel('x','FontSize',18);
        ylabel('h','FontSize',18);
        set(gca, 'FontSize', 18);
        hold off;
        
        subplot(2,2,2); %U Solution
        umax=max(max(max(u_movie)));
        umin=min(min(min(u_movie)));
        for e=1:nelem
            for j=1:ngl
                x(j)=coord(j,e);
                y(j)=u_movie(j,e,i);
            end
            plot_handle=plot(x,y,'r-');
            set(plot_handle,'LineWidth',2);
            hold on;
        end
        axis([xmin xmax umin umax]);
        xlabel('x','FontSize',18);
        ylabel('U','FontSize',18);
        set(gca, 'FontSize', 18);
        hold on;
        
        %Plot Elements
        for e=1:nelem
            x1=coord(1,e);
            x2=coord(ngl,e);
            y1=u_movie(1,e,i);
            y2=u_movie(ngl,e,i);
            plot(x1,y1,'ko');
            plot(x2,y2,'ko');
        end
        hold off;
        
        subplot(2,2,3); %L2 Norms
        ymax=max(max(L2_norm));
        ymin=min(min(L2_norm));
        xmax=max(time_movie);
        xmin=min(time_movie);
        xt=time_movie(1:i);
        h_norm=L2_norm(1,1:i);
        plot_handle=plot(xt,h_norm,'b-','LineWidth',2);
        u_norm=L2_norm(2,1:i);
        hold on;
        plot_handle=plot(xt,u_norm,'k-','LineWidth',2);
        xlabel('Time','FontSize',18);
        ylabel('L^2 Norm','FontSize',18);
        set(gca, 'FontSize', 18);
        legend('||h||_2','||u||_2');
        axis([ xmin xmax ymin ymax]);
        hold on;
        
        subplot(2,2,4); %Mass Conservation
        ymax=max(mass_movie);
        ymin=min(mass_movie);
        xmax=max(time_movie);
        xmin=min(time_movie);
        xt=time_movie(1:i);
        yt=mass_movie(1:i);
        plot_handle=plot(xt,yt,'r-','LineWidth',2);
        xlabel('Time','FontSize',18);
        ylabel('\Delta M','FontSize',18);
        set(gca, 'FontSize', 18);
        axis([ xmin xmax ymin ymax]);
        hold on;
%         xt=time_movie(1:i);
%         yt=energy_movie(1:i);
%         plot_handle=plot(xt,yt,'b-','LineWidth',2);
%         legend('\Delta h','\Delta U');
%         hold on;
        
        M_i=getframe(gcf);
        M(i)=M_i;
        pause(0.01);
        hold off;
    end
    if store_movie == 1
        file_movie=[main_text '_e' num2str(nelem) '_p' num2str(nop) '_q' num2str(noq) '_diss' num2str(diss) '.avi'];
        movie2avi(M,file_movie,'fps',5);
    end

elseif (iplot_movie == 2)
    
    figure; %error norms
    for i=1:iframe
        ymax=max(max(L2_norm));
        ymin=min(min(L2_norm));
        xmax=max(time_movie);
        xmin=min(time_movie);
        xt=time_movie(1:i);
        h_norm=L2_norm(1,1:i);
        plot_handle=plot(xt,h_norm,'b-','LineWidth',2);
        u_norm=L2_norm(2,1:i);
        hold on;
        plot_handle=plot(xt,u_norm,'k--','LineWidth',2);
        xlabel('Time','FontSize',18);
        ylabel('L^2 Norm','FontSize',18);
        set(gca, 'FontSize', 18);
        legend('||h||_2','||u||_2');
        axis([ xmin xmax ymin ymax]);
        hold on;
        M_i=getframe(gcf);
%        M1(i)=M_i;
        pause(0.01);
        hold off;
    end
    
    figure; %Mass loss
    for i=1:iframe
        ymax=max(mass_movie);
        ymin=min(mass_movie);
        xmax=max(time_movie);
        xmin=min(time_movie);
        xt=time_movie(1:i);
        yt=mass_movie(1:i);
        plot_handle=plot(xt,yt,'r-','LineWidth',2);
        xlabel('Time','FontSize',18);
        ylabel('\Delta M','FontSize',18);
        set(gca, 'FontSize', 18);
        axis([ xmin xmax ymin ymax]);
        hold on;
        M_i=getframe(gcf);
%        M(i)=M_i;
        pause(0.01);
        hold off;
    end
end
