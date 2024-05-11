%---------------------------------------------------------------------%
%---------------------------------------------------------------------%
%This code computes the 1D Shallow Water Equations using CG with either 
%LGL-only  or LGL+LGR points.
%It solves either linear or nonlinear SWE using delta_nl=0,1.
%Also, it can solve advection by using icase==0. All other cases uses SWE.
%The time-integration is accomplished via 2nd, 3rd Order, or 4th order RK
%Written by F.X. Giraldo on 12/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
%---------------------------------------------------------------------%
%Notes: LGR is better with increasing N_LGR. For LGL, LGR_BC=1 makes it 
%better than LGR. With LGR_BC=0 LGR is better than LGL (5x-10x). It matters 
%what the exponent and the type of sponge are. E.g. tanh behaves better than sine.
%The results are consistent for both linear and nonlinear.

%For LGL with LGR_BC=0, we do better with nrbc_right=2 => outflow BC. No
%Sponge required!
%---------------------------------------------------------------------%

clear all; 
close all;

tic

addpath('fxg') %New LGL and LGR points and bases

%------------------DO NOT CHANGE-------------------%
diss=1; %=1 dissipation (Rusanov Flux) and =0 no dissipation (Central Flux)
stages=4; %RK1, RK2, RK3-SSP, or RK4
Courant_max=0.25; %This is what controls the time-step
form_method='strong'; %weak or strong

%sponge variables
sponge_amp=1; %factor that multiplies sponge. >1 is stronger
sponge_exponent=4.0;
sponge_shape=1; %1=sin, 2=tanh, 3=shifted hyperbola (Modave)
sponge_method=1; %1=modify solution (hard reset) as in Giraldo SISC 2013, 2=gently via RHS operator in CREATE_RHS_VOLUME

%viscosity variables
lvisc=1;  
tau=0e1; %the larger the better for CG.
visc=1e-2;
C1=-0.25;
C2=0.25;

%Computes absolute maxs throughtout simulation
store_qmax=1;
%------------------DO NOT CHANGE-------------------%

%----------------------------------Plotting Variables----------------------%
%Plotting variables
nplots=10; %plotting variable - How many Frames
iprint=1000; %how often to print time-step info
plot_movie=1;
plot_figures=1;
plot_figures_title=0;
plot_elements=1;
plot_qmax=1;
plot_visc=0;
plot_mass=0;
plot_bathymetry=0;
plot_partial_domain=0; %2=only plot right-half of domain, 1=only plot part of domain, 0=plot entire domain
plot_reference_solution=1;
%----------------------------------Plotting Variables----------------------%

%------------------Allowed to CHANGE-------------------%
%Defines Test Case
time_final=2; %final time in revolutions: use =1 for NFBC and =10 for NRBC
delta_nl=0; %=0 linear and =1 nonlinear
icase=1; %=0 is a Gaussian exp(-16r^2) Advection as used in many other sample problems in Giraldo 2020 book
         %=1 is a Gaussian exp(-16r^2) as used in many other sample problems in Giraldo 2020 book
         %=2 is a Gaussian exp(-16r^2) but amp=0.1 and hmean=1.0
experiment_method=3; %0=LGL with NFBC
                     %1=Reference Solution: LGL with NRBC Standard: portion of the domain is sponge with sponge_amp=0
                     %2=LGL with NRBC Standard: portion of the domain is sponge (same polynomial degree)
                     %3=LGR NRBC 
                     %4=LGL NRBC: similar to LGR in that the sponge zone is scaled by a factor NOP_LGL 
                     %5=LGL NRBC: last element scaled to be same size as LGR element
                     %6=LGL NRBC: last element scaled to be same size as LGR element with homogeneous Dirichlet BC
                     %7=LGL NRBC: with outflow/extrapolation BC and No Sponge => Works really well!
                     %8=LGL NRBC: with outflow/extrapolation BC and No Sponge and last element scaled to be same size as LGR element
%------------------Allowed to CHANGE-------------------%

%----------------------------------Experiment Definitions----------------------%
%----LGL NFBC configuration
if (experiment_method==0)
    nelem_LGL=32; %Number of Elements
    nelem_LGR=1;
    nop_LGL=4;    %Polynomial Order
    nop_LGR=4;
    nelem=nelem_LGL+nelem_LGR;
    LGR_artificial_layer=0; %=1 increases size of element =0 keeps the same max/mins of domain
    LGR_scale=1; %1=scale and 0=do not scale
    LGR_scale_factor=1; %factor to scale from xgr_scale
    LGR_basis='LGL'; %LGL or LGR
    LGR_BC=0; %1=homogeneous Dirichlet BC, 0=turn it off
    lsponge=0;
    xmax_plot_type=1; %1=xmax_LGR, 2=xmax_LGL
    nrbc_right=1; %1= no-flux/hard-wall and 2=nrbc/extrapolation

%----LGL NRBC: Standard Method
elseif (experiment_method==1)
    nelem_LGL=32; %Number of Elements
    nelem_LGR=8;
    nop_LGL=4;    %Polynomial Order
    nop_LGR=4;
    nelem=nelem_LGL+nelem_LGR;
    LGR_artificial_layer=1; %=1 increases domain size to be ; =0 keeps the same max/mins of domain
    LGR_scale=1; %1=scale and 0=do not scale
    LGR_scale_factor=1; %factor to scale from xgr_scale
    LGR_basis='LGL'; %LGL or LGR
    LGR_BC=0; %1=homogeneous Dirichlet BC, 0=turn it off
    lsponge=1;
    sponge_amp=0; 
    xmax_plot_type=1; %1=xmax_LGR, 2=xmax_LGL
    nrbc_right=1; %1= no-flux/hard-wall and 2=nrbc/extrapolation

%----LGL NRBC: Standard Method
elseif (experiment_method==2)
    nelem_LGL=32; %Number of Elements
    nelem_LGR=8;
    nop_LGL=4;    %Polynomial Order
    nop_LGR=4;
    nelem=nelem_LGL+nelem_LGR;
    LGR_artificial_layer=1; %=1 increases domain size to be ; =0 keeps the same max/mins of domain
    LGR_scale=1; %1=scale and 0=do not scale
    LGR_scale_factor=1; %factor to scale from xgr_scale
    LGR_basis='LGL'; %LGL or LGR
    LGR_BC=0; %1=homogeneous Dirichlet BC, 0=turn it off
    lsponge=1;
    xmax_plot_type=1; %1=xmax_LGR, 2=xmax_LGL
    nrbc_right=1; %1= no-flux/hard-wall and 2=nrbc/extrapolation

%----LGR NRBC: 
elseif (experiment_method==3)
    nelem_LGL=32; %Number of Elements
    nelem_LGR=1;
    nop_LGL=4;    %Polynomial Order
    nop_LGR=32;
    nelem=nelem_LGL+nelem_LGR;
    LGR_artificial_layer=1; %=1 increases size of domain; =0 keeps the same max/mins of domain
    LGR_scale=1; %1=scale and 0=do not scale
    LGR_scale_factor=8; %factor to scale from xgr_scale
    LGR_basis='LGR'; %LGL or LGR
    LGR_BC=1; %1=homogeneous Dirichlet BC, 0=turn it off
    lsponge=1;
    xmax_plot_type=1; %1=xmax_LGR, 2=xmax_LGL
    nrbc_right=1; %1= no-flux/hard-wall and 2=nrbc/extrapolation

%----LGL NRBC: Similar to LGR in that last element is scaled by NOP_LGR
elseif (experiment_method==4)
    nelem_LGL=32; %Number of Elements
    nelem_LGR=1;
    nop_LGL=4;    %Polynomial Order
    nop_LGR=32;
    nelem=nelem_LGL+nelem_LGR;
    LGR_artificial_layer=1; %=1 increases size of element =0 keeps the same max/mins of domain
    LGR_scale=1; %1=scale and 0=do not scale
    LGR_scale_factor=8;%nop_LGR; %factor to scale from xgr_scale
    LGR_basis='LGL'; %LGL or LGR
    LGR_BC=0; %1=homogeneous Dirichlet BC, 0=turn it off
    lsponge=1;
    xmax_plot_type=1; %1=xmax_LGR, 2=xmax_LGL
    nrbc_right=1; %1= no-flux/hard-wall and 2=nrbc/extrapolation

%----LGL NRBC: Similar to LGR in that last element is scaled similarly to LGR element size 
elseif (experiment_method==5)
    nelem_LGL=32; %Number of Elements
    nelem_LGR=1;
    nop_LGL=4;    %Polynomial Order
    nop_LGR=4;
    nelem=nelem_LGL+nelem_LGR;
    LGR_artificial_layer=1; %=1 increases size of element =0 keeps the same max/mins of domain
    LGR_scale=1; %1=scale and 0=do not scale
    LGR_scale_factor=1e3; %factor to scale from xgr_scale
    LGR_basis='LGL'; %LGL or LGR
    LGR_BC=0; %1=homogeneous Dirichlet BC, 0=turn it off
    lsponge=1;
    xmax_plot_type=2; %1=xmax_LGR, 2=xmax_LGL
    nrbc_right=1; %1= no-flux/hard-wall and 2=nrbc/extrapolation

%----LGL NRBC: Similar to LGR in that last element is scaled similarly to LGR element size with Homogeneous Dirichlet BC
elseif (experiment_method==6)
    nelem_LGL=32; %Number of Elements
    nelem_LGR=1;
    nop_LGL=4;    %Polynomial Order
    nop_LGR=4;
    nelem=nelem_LGL+nelem_LGR;
    LGR_artificial_layer=1; %=1 increases size of element =0 keeps the same max/mins of domain
    LGR_scale=1; %1=scale and 0=do not scale
    LGR_scale_factor=1e4; %factor to scale from xgr_scale
    LGR_basis='LGL'; %LGL or LGR
    LGR_BC=1; %1=homogeneous Dirichlet BC, 0=turn it off
    lsponge=1;
    xmax_plot_type=2; %1=xmax_LGR, 2=xmax_LGL
    nrbc_right=1; %1= no-flux/hard-wall and 2=nrbc/extrapolation

%----LGL NRBC: with Outflow BC and NO SPONGE
elseif (experiment_method==7)
    nelem_LGL=32; %Number of Elements
    nelem_LGR=1;
    nop_LGL=4;    %Polynomial Order
    nop_LGR=4;
    nelem=nelem_LGL+nelem_LGR;
    LGR_artificial_layer=1; %=1 increases size of element =0 keeps the same max/mins of domain
    LGR_scale=1; %1=scale and 0=do not scale
    LGR_scale_factor=1; %factor to scale from xgr_scale
    LGR_basis='LGL'; %LGL or LGR
    LGR_BC=0; %1=homogeneous Dirichlet BC, 0=turn it off
    lsponge=1;
    sponge_amp=0;
    xmax_plot_type=1; %1=xmax_LGR, 2=xmax_LGL
    nrbc_right=2; %1= no-flux/hard-wall and 2=nrbc/extrapolation

%----LGL NRBC: Similar to LGR with Outflow BC and NO SPONGE and last element scaled similarly to LGR element size
elseif (experiment_method==8)
    nelem_LGL=32; %Number of Elements
    nelem_LGR=1;
    nop_LGL=4;    %Polynomial Order
    nop_LGR=4;
    nelem=nelem_LGL+nelem_LGR;
    LGR_artificial_layer=1; %=1 increases size of element =0 keeps the same max/mins of domain
    LGR_scale=1; %1=scale and 0=do not scale
    LGR_scale_factor=1e4; %factor to scale from xgr_scale
    LGR_basis='LGL'; %LGL or LGR
    LGR_BC=0; %1=homogeneous Dirichlet BC, 0=turn it off
    lsponge=1;
    sponge_amp=0;
    xmax_plot_type=2; %1=xmax_LGR, 2=xmax_LGL
    nrbc_right=2; %1= no-flux/hard-wall and 2=nrbc/extrapolation
end
%----------------------------------Experiment Definitions----------------------%

%For standing wave, needs to be linear and with no viscosity
%time_final=0.5 is when u will be perfectly zero.
if (icase == 4)
    lvisc=0;
    delta_nl=0;
end

%Store Constants
min_visc_elem=0; max_visc_elem=0;

%Record which diffusion is being used
if (lvisc == 0)
   diffusion = 'NoDiffusion';
   visc=0;
else
    if (C1 > 0 && C2 > 0)
        diffusion = 'SGS';
    else
        diffusion = 'Constant';
    end
end
%main_text=diffusion;

%Store Constants
ngl_LGL=nop_LGL + 1;
ngl_LGR=nop_LGR + 1;
ngl=zeros(nelem,1);
for e=1:nelem_LGL
    ngl(e)=ngl_LGL;
end
for e=nelem_LGL+1:nelem
    ngl(e)=ngl_LGR;
end

%----------------------------------------------------------------------------------------------%
%---------------------Part A: Students add these Functions---------------------------------------------%
%Compute Interpolation/Integration Points and Basis Functions
[xgl,wgl]=legendre_gauss_lobatto(ngl_LGL);
[psil,dpsil] = lagrange_basis3(ngl_LGL,ngl_LGL,xgl,xgl);
if strcmp(LGR_basis,'LGL')
    [xgr,wgr]=legendre_gauss_lobatto(ngl_LGR);
    [psir,dpsir] = lagrange_basis3(ngl_LGR,ngl_LGR,xgr,xgr);
elseif strcmp(LGR_basis,'LGR')
    ipoints=5; %LGR
    [xgr,wgr]=laguerre_gauss_radau_eigenvalue_fxg(ngl_LGR);
    [psir,dpsir,~,~] = lagrange_basis_laguerre_v3(ngl_LGR,ngl_LGR,xgr); 
end
main_text=LGR_basis;

%Store Basis Function derivatives in one element-wise array
dpsi=unite_dpsi(dpsil,dpsir,ngl,nelem,nelem_LGL);
%---------------------Part A: Students add these Functions---------------------------------------------%
%----------------------------------------------------------------------------------------------%

%Create Grid
[coord,coord_cg,intma,jac,wnq,npoin,xmin_grid,xmax_grid] = create_grid(ngl,ngl_LGR,nelem,nelem_LGL,xgl,xgr,wgl,wgr,icase,LGR_basis,LGR_scale,LGR_scale_factor,LGR_artificial_layer,lsponge,delta_nl);
dx_min=min(coord(2,1)-coord(1,1),coord(2,nelem)-coord(1,nelem)); %min of 1st and last elements

%Compute Exact Solution
time=0;
[qe,qb,gravity] = exact_solution(intma,coord,npoin,nelem,ngl,time,icase);

%Estimate time-step
dt=compute_dt(qe,qb,intma,ngl,nelem,icase,Courant_max,dx_min,gravity);

%Compute Initial Mass
[mass0,energy0]=compute_mass(qe,qb,intma,ngl,nelem,wnq,jac);

%----------------------------------------------------------------------------------------------%
%---------------------Part A: Students add these Functions---------------------------------------------%
%Create Local/Element Mass
Mmatrix = create_Mmatrix(intma,jac,npoin,nelem,ngl,wnq);
%---------------------Part A: Students add these Functions---------------------------------------------%
%----------------------------------------------------------------------------------------------%

%Klemp-Lily sponge
[Igamma,Igamma_volume,xmax_LGL,xt,xd]=create_sponge(qb,gravity,coord,coord_cg,npoin,ngl_LGL,nelem_LGL,lsponge,sponge_method,sponge_shape,sponge_amp,sponge_exponent);

%Define Plotting Maxima
if (xmax_plot_type==1)
    xmax_plot=xt;
elseif (xmax_plot_type==2)
    xmin=min(coord_cg);
    dx_plot=(xd-xmin)/nelem_LGL;
    xmax_plot=xmax_LGL+dx_plot;
end

%Initialize State Vector
qp=qe; q0=qe; q1=qe;
rhs=zeros(npoin,2);
rhs_rk4=zeros(npoin,2,4);
visc_elem=zeros(npoin,2);
iframe=0;

%Initialize Movie arrays
ntime=time_final/dt;
iplot=round(ntime/nplots);
ivector=0;
nvector=1000000;
h_movie=zeros(npoin,2*nplots);
u_movie=zeros(npoin,2*nplots);
time_movie=zeros(2*nplots,1);
mass_movie=zeros(2*nplots,1);
energy_movie=zeros(2*nplots,1);
time_vector=zeros(nvector,1);
qmax_vector=zeros(nvector,2);
mass_vector=zeros(nvector,1);

%Compute RK Time-Integration Coefficients
[a0,a1,beta] = compute_ti_coefficients(stages);

%-----------------------------------Time Integration
itime=0;
while time < time_final
   %Modify DT
   [dt,u_max]=compute_dt(qe,qb,intma,ngl,nelem,icase,Courant_max,dx_min,gravity);
   courant=u_max*dt/dx_min;
   time=time + dt;
   itime=itime + 1;

   %Make sure we finish at Time_Final
   if time>time_final
       time=time - dt;
       dt=time_final - time;
       time=time + dt;
   end
   
   for s=1:stages
      %------------------------------------Inviscid Operators------------------------------------------------%
      %---------------------Part B: Students add these Functions---------------------------------------------%
      %Create RHS Matrix
      rhs = create_rhs_volume(qp,qb,intma,jac,npoin,nelem,ngl,wnq,dpsi,gravity,delta_nl,Igamma_volume,icase,form_method); %Rayleigh damping in RHS
      %rhs = create_rhs_volume_weak(qp,qb,intma,jac,npoin,nelem,ngl,wnq,dpsi,gravity,delta_nl,Igamma_volume,icase);
      % rhs = create_rhs_volume_v2(qp,qb,intma,jac,npoin,nelem,ngl,wnq,dpsi,gravity,delta_nl,Igamma_volume,icase); %Rayleigh damping in RHS
      % max(rhs(:,1))
      % rhs
      if (icase >= 1)
        rhs = create_rhs_flux(rhs,qp,qb,intma,nelem,ngl,diss,gravity,delta_nl,lsponge,nrbc_right,form_method);
      end
      %Apply Communicator
      rhs = create_rhs_communicator(rhs,Mmatrix,npoin);
      %---------------------Part B: Students add these Functions---------------------------------------------%
      %------------------------------------Inviscid Operators------------------------------------------------%

      %------------------------------------Viscous Operators-------------------------------------------------%
      %---------------------Part C: Students add these Functions---------------------------------------------%
      %Local Adaptive Viscosity
      visc_elem(:,1)=visc;
      visc_elem(:,2)=visc;
      if (lvisc == 1)
        [visc_elem,min_visc_elem,max_visc_elem] = compute_LAV_viscosity(rhs,qp,qb,intma,coord,npoin,nelem,ngl,gravity,visc,C1,C2);
        q_visc=qp;
        rhs_visc = create_Lmatrix_IBP_General(intma,jac,npoin,nelem,ngl,wnq,dpsi,q_visc,visc_elem);
        rhs_visc = create_Flux_SIPDG_General(rhs_visc,intma,coord,nelem,ngl,dpsi,q_visc,tau,visc_elem,lsponge);
        rhs_visc = create_rhs_communicator(rhs_visc,Mmatrix,npoin);
        q_visc=rhs_visc;
        rhs=rhs + rhs_visc;
      end
      %---------------------Part C: Students add these Functions---------------------------------------------%
      %------------------------------------Viscous Operators-------------------------------------------------%

      %New Stage Value
      for I=1:npoin
          qp(I,:)=( a0(s)*q0(I,:) + a1(s)*q1(I,:) + dt*beta(s)*rhs(I,:) )*Igamma(I);
      end
      
      %Homogeneous Dirichlet BC
      if (LGR_BC==1)
          qp(npoin,:)=0; 
      end

      %Update
      rhs_rk4(:,:,s)=rhs(:,:);
      if (icase == 0)
          qp(1,1)=0; %homogeneous dirichlet BC on the left
          qp(:,2)=q0(:,2);
      end
      q1=qp;

   end %s
   
   %RK Update Solution
   if (stages == 4)
       for I=1:npoin
           qp(I,1)=(q0(I,1) + dt/6.0*( rhs_rk4(I,1,1) + 2*rhs_rk4(I,1,2)...
                   + 2*rhs_rk4(I,1,3) + rhs_rk4(I,1,4) ))*Igamma(I);
           qp(I,2)=(q0(I,2) + dt/6.0*( rhs_rk4(I,2,1) + 2*rhs_rk4(I,2,2)...
                   + 2*rhs_rk4(I,2,3) + rhs_rk4(I,2,4) ))*Igamma(I);
       end
   end
         
   %Homogeneous Dirichlet BC
   if (LGR_BC==1)
       qp(npoin,:)=0; 
   end

   %Update Q
   if (icase == 0)
      qp(1,1)=0;
      qp(:,2)=q0(:,2);
   end
   q0=qp;
   
   %Store Movie Frames
   if (plot_movie == 1 && mod(itime,iplot) == 0)
      [h_movie,u_movie,time_movie,mass_movie,energy_movie,iframe] = store_movie_frames(qp,qb,intma,jac,wnq,nelem,ngl,h_movie,u_movie,time_movie,mass_movie,energy_movie,iframe,mass0,energy0,time);           
   end 
   
   %Store Max values
   if (store_qmax == 1)
       ivector=ivector+1;
       time_vector(ivector)=time;
       if (plot_partial_domain==0)
            qmax_vector(ivector,1)=max(abs(qp(:,1)));
            qmax_vector(ivector,2)=max(abs(qp(:,2)));
        elseif (plot_partial_domain==1)
            xmax=-xmin_grid;
            xmin=xmin_grid;
            qmax_vector(ivector,1)=max(abs(qp(:,1)));
            qmax_vector(ivector,2)=max(abs(qp(:,2)));
        elseif (plot_partial_domain==2)
            hmax=-1000;
            umax=-1000;
            for I=1:npoin
                x=coord_cg(I);
                if (x>0) 
                    hmax=max(hmax,abs(qp(I,1)));
                    umax=max(umax,abs(qp(I,2)));
                end
            end
            qmax_vector(ivector,1)=hmax;
            qmax_vector(ivector,2)=umax;
        end
        [mass,energy]=compute_mass(qp,qb,intma,ngl,nelem,wnq,jac);
        mass_vector(ivector)=abs(mass-mass0);%/mass0;
   end 

   %Print to Screen
   if mod(itime,iprint) == 0
       % disp(['itime =  ',num2str(itime),' time = ', num2str(time),' dt = ', num2str(dt),' courant = ', num2str(courant),' |h|_max = ', num2str(max(abs(qp(:,1)))),' |u|_max = ', num2str(max(abs(qp(:,2))))]);
        disp(['itime =  ',num2str(itime),' time = ', num2str(time),' dt = ', num2str(dt),' courant = ', num2str(courant),' |h|_max = ', num2str(qmax_vector(ivector,1)),' |u|_max = ', num2str(qmax_vector(ivector,2))]);

   end

end %itime
%----------------------------------End Time Integration


%----------------------------------Begin Plotting Results
if (plot_movie == 1)
    
    figure;
    for i=1:iframe
        subplot(3,1,1); %---------------------------------------H Solution
        hmax=max(h_movie(:));
        hmin=min(h_movie(:));
        if plot_bathymetry==1
            hmin=min(-qb);
        end
        if (plot_partial_domain==0)
            xmax=xmax_plot;
            xmin=min(coord_cg);
        elseif (plot_partial_domain==1)
            xmax=-xmin_grid;
            xmin=xmin_grid;
        elseif (plot_partial_domain==2)
            xmax=xmax_plot;
            xmin=0;
        end
        plot(coord_cg,h_movie(:,i),'r-','LineWidth',2);
        hold on;
        if plot_bathymetry==1
            plot(coord_cg,-qb,'b-','LineWidth',2);
        end
        axis([xmin xmax hmin hmax]);
        title_text=[main_text ' visc = ' num2str(visc) ', Ne = ' num2str(nelem) ', NLGL = ' num2str(nop_LGL) ', NLGR = ' num2str(nop_LGR) ', Time = ' num2str(time_movie(i))];
        title([title_text],'FontSize',18);
        hold on;
         
        %Plot Elements
        if (plot_elements == 1)
            for e=1:nelem
                i1=intma(1,e);
                i2=intma(ngl(e),e);
                x1=coord(1,e);
                x2=coord(ngl(e),e);
                y1=h_movie(i1,i);
                y2=h_movie(i2,i);
                plot(x1,y1,'ko');
                plot(x2,y2,'ko');
            end
        end

        %Plot LGL/LGR boundary
        e=nelem_LGL;
        I=intma(ngl(e),e);
        xx(1)=coord_cg(I);
        xx(2)=coord_cg(I);
        yy(1)=hmin;
        yy(2)=hmax;
        if (lsponge==1)
            plot(xx,yy,'b:','LineWidth',2.0);
        end

        %Labels and Title
        title([title_text],'FontSize',18);
        xlabel('x','FontSize',18);
        ylabel('h','FontSize',18);
        set(gca, 'FontSize', 18);
        hold off;
        
        subplot(3,1,2); %----------------------------------------U Solution
        umax=max(u_movie(:));
        umin=min(u_movie(:));
        plot(coord_cg,u_movie(:,i),'r-','LineWidth',2);
        axis([xmin xmax umin umax]);
        xlabel('x','FontSize',18);
        ylabel('U','FontSize',18);
        set(gca, 'FontSize', 18);
        hold on;
        
        %Plot Elements
        if (plot_elements == 1)
            for e=1:nelem
                i1=intma(1,e);
                i2=intma(ngl(e),e);
                x1=coord(1,e);
                x2=coord(ngl(e),e);
                y1=u_movie(i1,i);
                y2=u_movie(i2,i);
                plot(x1,y1,'ko');
                plot(x2,y2,'ko');
            end
        end
        
        %Plot LGL/LGR boundary
        e=nelem_LGL;
        I=intma(ngl(e),e);
        xx(1)=coord_cg(I);
        xx(2)=coord_cg(I);
        yy(1)=umin;
        yy(2)=umax;
        if (lsponge==1)
            plot(xx,yy,'b:','LineWidth',2.0);
        end

        %Labels and Title
        title([title_text],'FontSize',18);
        xlabel('x','FontSize',18);
        ylabel('U','FontSize',18);
        set(gca, 'FontSize', 18);
        hold off;

        subplot(3,1,3); %--------------------------------Mass Conservation
        ymax=max(mass_movie);
        %ymax=1e-15;
        ymin=min(mass_movie);
        xmax=max(time_movie);
        xmin=min(time_movie);
        xt=time_movie(1:i);
        yt=mass_movie(1:i);
        semilogy(xt,yt,'r-','LineWidth',2);
        xlabel('Time','FontSize',18);
        ylabel('\Delta M','FontSize',18);
        set(gca, 'FontSize', 18);
        axis([ xmin xmax ymin ymax]);
        hold on;
        
        M_i=getframe(gcf);
        M(i)=M_i;
        pause(0.01);
        hold off;
    end
end

if (plot_figures == 1)
        
    figure; %-------------------------------------------------------H Plot
    hmax=max(qp(:,1));
    hmin=min(qp(:,1));
    if plot_bathymetry==1
        hmin=min(-qb);
    end
    
    if (plot_partial_domain==0)
        xmax=xmax_plot;
        xmin=min(coord_cg);
    elseif (plot_partial_domain==1)
        xmax=-xmin_grid;
        xmin=xmin_grid;
    elseif (plot_partial_domain==2)
        xmax=xmax_plot;
        xmin=0;
        hmax=-1000;
        hmin=+1000;
        for i=1:npoin
            x=coord_cg(i);
            if (x>0) 
                hmax=max(hmax,qp(i,1));
                hmin=min(hmin,qp(i,1));
            end
        end
    end
    plot(coord_cg,qp(:,1),'r-','LineWidth',2);
    hold on;
    if plot_bathymetry==1
        plot(coord_cg,-qb,'b-','LineWidth',2);
    end
    if (plot_reference_solution==1)
        load('1D_SWE_Linear_Exact_t=2.mat','q_ref','coord_ref');
        plot(coord_ref,q_ref(:,1),'k:','LineWidth',2);
    end
    axis([xmin xmax hmin hmax]);
    title_text=[main_text ' visc = ' num2str(visc) ', Ne = ' num2str(nelem) ', NLGL = ' num2str(nop_LGL) ', NLGR = ' num2str(nop_LGR) ', Time = ' num2str(time)];
    hold on;

    %Plot Elements
    if (plot_elements == 1)
        for e=1:nelem
            i1=intma(1,e);
            i2=intma(ngl(e),e);
            x1=coord(1,e);
            x2=coord(ngl(e),e);
            y1=qp(i1,1);
            y2=qp(i2,1);
            plot(x1,y1,'ko');
            plot(x2,y2,'ko');
        end
    end

    %Plot LGL/LGR boundary
    e=nelem_LGL;
    I=intma(ngl(e),e);
    xx(1)=coord_cg(I);
    xx(2)=coord_cg(I);
    yy(1)=hmin;
    yy(2)=hmax;
    if (lsponge==1)
       plot(xx,yy,'b:','LineWidth',2.0);
    end

    %Labels and Title
    if plot_figures_title == 1
        title(title_text,'FontSize',18);
    end
    xlabel('x','FontSize',18);
    ylabel('h','FontSize',18);
    if plot_bathymetry==0
        ylabel('h_s','FontSize',18);
    end
    set(gca, 'FontSize', 18);
    hold off;
    
    figure; %-------------------------------------------------------U Plot
    hmax=max(qp(:,2));
    hmin=min(qp(:,2));
    if (icase == 0)
        hmin=0;
    end
    if (plot_partial_domain==0)
            xmax=xmax_plot;
            xmin=min(coord_cg);
    elseif (plot_partial_domain==1)
        xmax=-xmin_grid;
        xmin=xmin_grid;
    elseif (plot_partial_domain==2)
        xmax=xmax_plot;
        xmin=0;
        hmax=-1000;
        hmin=+1000;
        for i=1:npoin
            x=coord_cg(i);
            if (x>0) 
                hmax=max(hmax,qp(i,2));
                hmin=min(hmin,qp(i,2));
            end
        end
    end
    plot(coord_cg,qp(:,2),'r-','LineWidth',2);
    axis([xmin xmax hmin hmax]);
    title_text=[main_text ' visc = ' num2str(visc) ', Ne = ' num2str(nelem) ', NLGL = ' num2str(nop_LGL) ', NLGR = ' num2str(nop_LGR) ', Time = ' num2str(time)];
    hold on;
    
    %Plot Elements
    if (plot_elements == 1)
        for e=1:nelem
            i1=intma(1,e);
            i2=intma(ngl(e),e);
            x1=coord(1,e);
            x2=coord(ngl(e),e);
            y1=qp(i1,2);
            y2=qp(i2,2);
            plot(x1,y1,'ko');
            plot(x2,y2,'ko');
        end
    end

    %Plot LGL/LGR boundary
    e=nelem_LGL;
    I=intma(ngl(e),e);
    xx(1)=coord_cg(I);
    xx(2)=coord_cg(I);
    yy(1)=hmin;
    yy(2)=hmax;
    if (lsponge==1)
       plot(xx,yy,'b:','LineWidth',2.0);
    end

    %Labels and Title
    if plot_figures_title == 1
        title(title_text,'FontSize',18);
    end
    xlabel('x','FontSize',18);
    ylabel('U','FontSize',18);
    set(gca, 'FontSize', 18);
    hold off;
end %plot figures

%Plot Viscosity
if (plot_visc == 1)
    figure; %-------------------------------------------------VISC Plot
    hmax=max(visc_elem(:,1));
    hmin=0;
    xmax=xmax_plot;
    xmin=min(coord_cg);
    plot(coord_cg,visc_elem(:,1),'r-','LineWidth',2);
    axis([xmin xmax hmin hmax]);
    title_text=[main_text ' visc = ' num2str(visc) ', Ne = ' num2str(nelem) ', NLGL = ' num2str(nop_LGL) ', NLGR = ' num2str(nop_LGR) ', Time = ' num2str(time)];
    hold on;

    %Plot Elements
    if (plot_elements == 1)
        for e=1:nelem
            i1=intma(1,e);
            i2=intma(ngl(e),e);
            x1=coord(1,e);
            x2=coord(ngl(e),e);
            y1=visc_elem(i1,1);
            y2=visc_elem(i2,1);
            plot(x1,y1,'ko');
            plot(x2,y2,'ko');
        end
    end

    if plot_figures_title == 1
        title(title_text,'FontSize',18);
    end
    xlabel('x','FontSize',18);
    ylabel('VISC','FontSize',18);
    set(gca, 'FontSize', 18);
    hold off;
end %lvisc

%Plot Maximum Values at the End
if (plot_qmax==1)
    tt=zeros(ivector,1);
    qq=zeros(ivector,2);
    for i=1:ivector
        tt(i)=time_vector(i);
        qq(i,1)=qmax_vector(i,1);
        qq(i,2)=qmax_vector(i,2);
    end
    figure;
    hmax=max(qq(:));
    hmin=min(qq(:));
    xmax=max(tt(:));
    xmin=min(tt(:));
    semilogy(tt,qq(:,1),'r-','LineWidth',2);
    hold on;
    semilogy(tt,qq(:,2),'b-','LineWidth',2);
    axis([xmin xmax hmin hmax]);
    %axis([1 10 1e-11 1e0]);
    % axis([1 10 hmin hmax]);
    title_text=[main_text,' Ne = ' num2str(nelem) ', NLGL = ' num2str(nop_LGL) ', NLGR = ' num2str(nop_LGR)];
    hold on;

    legend('h','U')
    if plot_figures_title == 1
        title(title_text,'FontSize',18);
    end
    xlabel('time','FontSize',18);
    ylabel('max(|q(x,t^n)|)','FontSize',18);
    set(gca, 'FontSize', 18);
end

%Plot Mass with Time
if (plot_mass==1)
    tt=zeros(ivector,1);
    qq=zeros(ivector,1);
    for i=1:ivector
        tt(i)=time_vector(i);
        qq(i)=mass_vector(i);
    end
    figure;
    hmax=max(qq);
    hmin=min(qq);
    xmax=max(tt);
    xmin=min(tt);
    semilogy(tt,qq,'r-','LineWidth',2);
    hold on;
    axis([xmin xmax hmin hmax]);
    title_text=[main_text,' Ne = ' num2str(nelem) ', NLGL = ' num2str(nop_LGL) ', NLGR = ' num2str(nop_LGR)];
    hold on;
    if plot_figures_title == 1
        title(title_text,'FontSize',18);
    end
    xlabel('Time','FontSize',18);
    ylabel('\Delta M','FontSize',18);
    set(gca, 'FontSize', 18);
end