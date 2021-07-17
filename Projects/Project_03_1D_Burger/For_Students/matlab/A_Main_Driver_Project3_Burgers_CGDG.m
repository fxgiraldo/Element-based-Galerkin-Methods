%---------------------------------------------------------------------%
%This code contains the Student template for solving the 1D Burgers Equations 
%using either CG or DG method with either LG or LGL points for interpolation and integration.
%
%The time-integration is accomplished via 2nd, 3rd Order, or 3rd Order 4-stage
%RK.
%
%The strategy followed is very similar to that described in Alg. 7.6 except 
%that we do not build the Flux matrix explicitly. This is all done in 
%CREATE_RHS_VOLUME and CREATE_RHS_FLUX.
%
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
%--------------------------ONLY CHANGE THESE--------------------------%
nelem=21; %Number of Elements
nop=6;    %Polynomial Order
time_final=0.4; %final time in revolutions
integration_points=1; %=1 for LGL and =2 for LG
integration_type=1; %=1 is inexact and =2 is exact
space_method='dg'; %cgc=CGc; cgd=CGd; dg=DG
%--------------------------ONLY CHANGE THESE--------------------------%
%For nop=2,4,6 works upto Time=10 (nop*nelem=120) dt=0.001. 

kstages=4; %RK1, RK2, RK3, or RK4
dt=0.001;
nplots=10; %plotting variable - Number of Frames ploted
iplot_movie=0;
store_movie=0;
iplot_figures=1;
iplot_elements=1;
iplot_modes=1;
%alpha=2.0/3.0; %1=conservative, 0=non-conservative
alpha=1; %1=conservative, 0=non-conservative
limit=0; %=0 no limiting, =1 Krivodonova Limiter, =2 Bound-Preserving Limiter
icase=1; %=1 is a Gaussian with flat bottom; 
xmu=0; %filtering strength: 1 is full strength and 0 is no filter
ifilter=0; %time-step frequency that the filter is applied.
filter_type=2; %=1 is Modal Hierarchical (no need for DSS)
               %and =2 is regular Legendre (better for DG)
diss=1; %=1 dissipation (Rusanov Flux) and =0 no dissipation (Central Flux)
% visc=6.6e-4;
visc=0e-4;
tau=0e3;
flux_method='rusanov'; %rusanov,roe, or energy
LDG_alpha=0.5; LDG_beta=1.0-LDG_alpha;

%Viscous Operators
if (strcmp(space_method,'dg')>0)
    elliptic_method='LDG';
else
    elliptic_method='SIP';
end

%Store Constants
eps=1e-15;
ntime=time_final/dt;
iplot=round(ntime/nplots);

%Store Constants
ngl=nop + 1;
method_text = [space_method];

%Compute i,e -> I pointer
I=0;
Ipointer=zeros(ngl,nelem);
for e=1:nelem
    for i=1:ngl
        I=I+1;
        Ipointer(i,e)=I;
    end
end

%Compute Interpolation and Integration Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);
if (integration_points == 1)
    integration_text=['LGL'];
    if (integration_type == 1)
        noq=nop;
    elseif (integration_type == 2)
        noq=3*nop/2;
    end
    nq=noq + 1;
    [xnq,wnq]=legendre_gauss_lobatto(nq);
elseif (integration_points == 2)
    integration_text=['LG'];
    noq=nop;
    nq=noq + 1;
    [xnq,wnq]=legendre_gauss(nq);
end
main_text=[space_method   ':' integration_text ];

%Compute Lagrange Polynomial and derivatives
[psi,dpsi] = lagrange_basis3(ngl,nq,xgl,xnq);
   
%Compute Filter Matrix
[f,vdm,vdm_inv] = filter_init(ngl,xgl,xmu,filter_type);

%Create Grid
[coord,intma_cg]=create_grid_dg(ngl,nelem,xgl,icase);
dx_min=coord(2,1)-coord(1,1);

%Form Global Matrix Pointer
npoin_cg=nop*nelem + 1;
npoin_dg=ngl*nelem;
intma=zeros(ngl,nelem);
if ( strcmp(space_method,'cgc') > 0 )
    npoin=npoin_cg;
    intma=intma_cg;
else %cgd or dg
    npoin=npoin_dg;
    ip=0;
    for e=1:nelem
        for i=1:ngl
            ip=ip+1;
            intma(i,e)=ip;
        end
    end
end
if ( strcmp(space_method,'dg') > 0 )
    npoin_cg=npoin;
    intma_cg=intma;
end

%Compute Exact Solution
time=0;
p_movie=zeros(npoin,nplots);
u_movie=zeros(npoin,nplots);
time_movie=zeros(nplots,1);
mass_movie=zeros(nplots,1);
energy_movie=zeros(nplots,1);
[qe] = exact_solution_dg(intma,coord,npoin,nelem,ngl,time,icase);
q_init=qe;

%Compute Initial Mass and Energy
[mass0,energy0] = compute_Mass_and_Energy(qe,intma,coord,wnq,psi,nelem,ngl,nq);

%Create Periodic BC Pointer
iperiodic=zeros(npoin,1);
for i=1:npoin
   iperiodic(i)=i;
end
%iperiodic(npoin)=iperiodic(1);

%Create Local/Element Mass
Mmatrix = create_Mmatrix(intma_cg,coord,npoin_cg,nelem,ngl,nq,wnq,psi);

%Initialize State Vector
qp=qe; q0=qe; q1=qe;
iframe=0;

%Initialize Arrays
rhs=zeros(npoin,2);
rhs_t=zeros(ngl,1);
rhs_rk4=zeros(npoin,2,4);
visc_elem=zeros(npoin,2);

%Compute RK Time-Integration Coefficients
[a0,a1,beta] = compute_ti_coefficients(kstages);

%Time Integration
for itime=1:ntime
   time=time + dt;
   u_max=-100000;
   for e=1:nelem
       for i=1:ngl
           ip=intma(i,e);
           u_max=max(u_max, abs(qp(ip,1)) );
       end
   end
   courant=u_max*dt/dx_min;
    
   %disp(['itime =  ',num2str(itime),' time = ', num2str(time),' courant = ', num2str(courant)]);
   for ik=1:kstages
      
      %Create RHS vector for Inviscid Operators
      %-----------------------------------------------------------------------------%
      %-------------------Student add your RHS functions here ----------------------%
      %rhs = create_rhs_volume(qp,intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,alpha);
      %rhs = create_rhs_flux(rhs,qp,intma,nelem,ngl,diss,flux_method);
      %-----------------------------------------------------------------------------%
      %-----------------------------------------------------------------------------%
      
      %Create RHS vector for Viscous Operators
      if (strcmp(elliptic_method,'SIP') > 0 )
          rhs = create_Lmatrix_SIPDG(rhs,intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,qp,tau,visc);
      elseif (strcmp(elliptic_method,'LDG') > 0 )
          rhs = create_Lmatrix_LDG(rhs,Mmatrix,intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,qp,LDG_alpha,LDG_beta,visc);      
      end
      
      %Apply Communicator
      %-----------------------------------------------------------------------------%
      %-------------------Student add your RHS functions here ----------------------%
      %rhs = create_rhs_communicator(rhs,space_method,intma,intma_cg,coord,Mmatrix,npoin,npoin_cg,nelem,ngl,nq,wnq,psi,dpsi,alpha,0);
      %-----------------------------------------------------------------------------%
      %-------------------Student add your RHS functions here ----------------------%
      
      %Solve System
      qp=a0(ik)*q0 + a1(ik)*q1 + dt*beta(ik)*rhs;
                  
      %Update
      rhs_rk4(:,:,ik)=rhs(:,:);
      q1=qp;
   end %ik
   
   %Do RK4 Solution
   if (kstages == 4)
       qp(:,:)=q0(:,:) + dt/6.0*( rhs_rk4(:,:,1) + 2*rhs_rk4(:,:,2)...
                      + 2*rhs_rk4(:,:,3) + rhs_rk4(:,:,4) );
   end
   
   %Filter Solution
   if (mod(itime,ifilter) == 0 && ifilter > 0)
      qp = apply_filter_dg(qp,f,intma,nelem,ngl);
   end
   
   %Limit Solution
   limit_element=zeros(nelem,1);
   if limit == 0 || limit == 1
       [limit_element,qp,qmodal] = limiter_nodal_v3(qp,vdm,vdm_inv,intma,nelem,ngl,limit);
   elseif limit == 2
       [limit_element,qp,qmodal] = limiter_nodal_bp(qp,q0,vdm,vdm_inv,intma,nelem,ngl);
   end
   
   %DSS Solution for CG after Limiting
   if (limit > 0 || ( mod(itime,ifilter) == 0 && ifilter > 0 ) )
       if ( strcmp(space_method,'cgd') > 0 || strcmp(space_method,'cgc') > 0) %CG
          qp = apply_dss(qp,intma,intma_cg,coord,Mmatrix,npoin,npoin_cg,nelem,ngl,nq,wnq,psi,dpsi,alpha,0);
       end
   end
   
   %Plot Modes
   if (limit >= 0 && iplot_modes == 1 && mod(itime,iplot) == 0)
       [test] = plot_modes(qmodal,Ipointer,limit_element,intma,coord,qp,nelem,ngl,nq,nop,noq,time,diss,main_text);           
   end

   %Update Q
   q0=qp;
   qb=q0;
   qb=0;
   if (iplot_movie == 1 && mod(itime,iplot) == 0)
      [p_movie,u_movie,time_movie,mass_movie,energy_movie,L2_norm,iframe] = store_movie_frames(qp,q0,intma,coord,wnq,psi,npoin,nelem,ngl,nq,p_movie,u_movie,time_movie,mass_movie,energy_movie,iframe,mass0,energy0,time,icase);           
   end 
   
end %itime

if (iplot_movie == 1)
    
    hmax=max(p_movie(:));
    hmin=min(p_movie(:));
    xmax=max(max(coord));
    xmin=min(min(coord));
    figure;
    for i=1:iframe
        subplot(3,1,1); %H Solution
%         hmax=max(p_movie(:,i));
%         hmin=min(p_movie(:,i));
%         xmax=max(max(coord));
%         xmin=min(min(coord));
        for e=1:nelem
            for j=1:ngl
                jp=intma(j,e);
                x(j)=coord(j,e);
                y(j)=p_movie(jp,i);
            end
            plot_handle=plot(x,y,'r-','LineWidth',2);
            hold on;
        end
        axis([xmin xmax hmin hmax]);
        title_text=[main_text ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time_movie(i))];
        title([title_text],'FontSize',18);
        hold on;
         
        %Plot Elements
        for e=1:nelem
            i1=intma(1,e);
            i2=intma(ngl,e);
            x1=coord(1,e);
            x2=coord(ngl,e);
            y1=p_movie(i1,i);
            y2=p_movie(i2,i);
            plot(x1,y1,'ko');
            plot(x2,y2,'ko');
        end

        title([title_text],'FontSize',18);
        xlabel('x','FontSize',18);
        ylabel('h','FontSize',18);
        set(gca, 'FontSize', 18);
        hold off;
                
        subplot(3,1,2); %L2 Norms
        ymax=max(max(L2_norm));
        ymin=min(min(L2_norm));
        tmax=max(time_movie);
        tmin=min(time_movie);
        tt=time_movie(1:i);
        h_norm=L2_norm(1,1:i);
        plot_handle=plot(tt,h_norm,'b-','LineWidth',2);
        u_norm=L2_norm(2,1:i);
        hold on;
        plot_handle=plot(tt,u_norm,'k-','LineWidth',2);
        xlabel('Time','FontSize',18);
        ylabel('L^2 Norm','FontSize',18);
        set(gca, 'FontSize', 18);
        legend('||h||_2','||u||_2');
        axis([ tmin tmax ymin ymax]);
        hold on;
        
        subplot(3,1,3); %Mass Conservation
        ymax=max(mass_movie);
        ymin=min(mass_movie);
        tmax=max(time_movie);
        tmin=min(time_movie);
        tt=time_movie(1:i);
        yt=mass_movie(1:i);
        plot_handle=plot(tt,yt,'r-','LineWidth',2);
        xlabel('Time','FontSize',18);
        ylabel('\Delta M','FontSize',18);
        set(gca, 'FontSize', 18);
%        axis([ tmin tmax ymin ymax]);
        hold on;
        tt=time_movie(1:i);
        yt=energy_movie(1:i);
        plot_handle=plot(tt,yt,'b-','LineWidth',2);
        legend('\Delta M','\Delta E');
        hold on;
        
        M_i=getframe(gcf);
        M(i)=M_i;
        pause(0.01);
        hold off;
    end
    if store_movie == 1
        file_movie=[main_text '_e' num2str(nelem) '_p' num2str(nop) '_q' num2str(noq) '_diss' num2str(diss) '.avi'];
        movie2avi(M,file_movie,'fps',5);
    end
end

if (iplot_figures)
    
    %Exact Solution
    [qe] = exact_solution_dg(intma,coord,npoin,nelem,ngl,time,icase);
    
    figure; %H Plot
    hmax=max(qp(:,1));
    hmin=min(qp(:,1));
    xmax=max(coord(:));
    xmin=min(coord(:));
    for e=1:nelem
        for i=1:ngl
            ip=intma(i,e);
            x(i)=coord(i,e);
            y(i)=qp(ip,1);
            ye(i)=q_init(ip,1);
        end
        plot_handle=plot(x,y,'r-','LineWidth',2);
        hold on;
        plot_handle=plot(x,ye,'k:','LineWidth',2);
    end
    axis([xmin xmax hmin hmax]);
    title_text=[main_text ' visc = ' num2str(visc) ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time)];
    title([title_text],'FontSize',18);
    hold on;

    %Plot Elements
    if (iplot_elements == 1)
        for e=1:nelem
            i1=intma(1,e);
            i2=intma(ngl,e);
            x1=coord(1,e);
            x2=coord(ngl,e);
            y1=qp(i1,1);
            y2=qp(i2,1);
            plot(x1,y1,'ko');
            plot(x2,y2,'ko');
        end
    end

    title([title_text],'FontSize',18);
    xlabel('x','FontSize',18);
    ylabel('u(x,t)','FontSize',18);
    set(gca, 'FontSize', 18);
    hold off;
            
end