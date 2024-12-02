%---------------------------------------------------------------------%
%What to show for nel=20, N=8, time_final=1.0, dt=0.001
%1) NES with no filter => unstable
%2) NES with a filter (xmu=0.05) => stable
%3) ES with no filter => stable
%4) NES with a filter (xmu=0.01) => unstable
%5) ES with a filter (xmu=0.01) => stable
%---------------------------------------------------------------------%
%This code computes the 1D Burgers Equation using either CG or 
%DG method with LGL points for interpolation and integration.
%The code can use Vanilla or Entropy-Stable DG.
%The time-integration is accomplished via 2nd, 3rd Order, or 3rd Order 4-stage
%RK.
%Written by F.X. Giraldo on 5/2021
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%------------------DO NOT CHANGE-------------------%
stages=4; %RK1, RK2, RK3, or RK4
Courant_max=0.5;
icase=1; %=1 Sine function; 
diss=1; %=1 dissipation (Rusanov Flux) and =0 no dissipation (Central Flux)
space_method='dg'; %cgc=CGc; cgd=CGd; dg=DG
cgdg_method='strong'; %weak or strong (strong/weak work for NES but only strong works for ES)
integration_points=1; %=1 for LGL and =2 for LG
integration_type=1; %=1 is inexact and =2 is exact
%------------------DO NOT CHANGE-------------------%

%----------------------------------Plotting Variables----------------------%
nplots=10; %plotting variable - Number of Frames ploted
iplot_movie=1;
store_movie=0;
iplot_figures=1;
iplot_elements=1;
iplot_modes=0;
%----------------------------------Plotting Variables----------------------%

%------------------Allowed to CHANGE-------------------%
time_final=1.0; %final time in revolutions
nelem=20; %Number of Elements
nop=8;    %Polynomial Order
experiment_method=2; %1=Entropy Stable flux with no filter and no limiter => works
                     %2=Standard flux with no filter and no limiter => doesn't work
                     %3=Standard flux with filter => works
                     %4=Standard flux with limiter => doesn't blow-up but wrong solution
%------------------Allowed to CHANGE-------------------%

%----------------------------------Experiment Definitions----------------------%
main_text='Experiment ?';
%----Entropy-Stable Flux => works
if (experiment_method==1)
    volume_flux_method='ES'; %ES = Entropy-Stable or NES = Not Entropy-Stable (Classical flux function)
    surface_flux_method='ES'; %ES = Entropy-Stable or NES = Not Entropy-Stable (Classical flux function)
    limit=0; %=0 no limiting, =1 Krivodonova Limiter, =2 Bound-Preserving Limiter
    xmu=0.05; %filtering strength: 1 is full strength and 0 is no filter
    ifilter=0; %time-step frequency that the filter is applied.
    filter_type=1; 
    main_text='Experiment 1';
%----Standard Flux => doesn't work
elseif (experiment_method==2)
    volume_flux_method='NES'; %ES = Entropy-Stable or NES = Not Entropy-Stable (Classical flux function)
    surface_flux_method='NES'; %ES = Entropy-Stable or NES = Not Entropy-Stable (Classical flux function)
    limit=0; %=0 no limiting, =1 Krivodonova Limiter, =2 Bound-Preserving Limiter
    xmu=0.05; %filtering strength: 1 is full strength and 0 is no filter
    ifilter=0; %time-step frequency that the filter is applied.
    filter_type=1; 
    main_text='Experiment 2';
%----Standard Flux with Filter => works
elseif (experiment_method==3)
    volume_flux_method='NES'; %ES = Entropy-Stable or NES = Not Entropy-Stable (Classical flux function)
    surface_flux_method='NES'; %ES = Entropy-Stable or NES = Not Entropy-Stable (Classical flux function)
    limit=0; %=0 no limiting, =1 Krivodonova Limiter, =2 Bound-Preserving Limiter
    xmu=1.0; %filtering strength: 1 is full strength and 0 is no filter
    ifilter=1; %time-step frequency that the filter is applied.
    filter_type=1; 
    main_text='Experiment 3';
%----Standard Flux with Limiter => works
elseif (experiment_method==4)
    volume_flux_method='NES'; %ES = Entropy-Stable or NES = Not Entropy-Stable (Classical flux function)
    surface_flux_method='NES'; %ES = Entropy-Stable or NES = Not Entropy-Stable (Classical flux function)
    limit=2; %=0 no limiting, =1 Krivodonova Limiter, =2 Bound-Preserving Limiter
    xmu=0.05; %filtering strength: 1 is full strength and 0 is no filter
    ifilter=0; %time-step frequency that the filter is applied.
    filter_type=1; 
    main_text='Experiment 4';
end

%Store Constants
eps=1e-15;
ivolume_integrate=0;
if (filter_type == 2)
    ivolume_integrate=1;
end

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
%main_text=[space_method  ':' cgdg_method ':' flux_method ];

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
if strcmp(space_method,'cgc')
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
if strcmp(space_method,'dg')
    npoin_cg=npoin;
    intma_cg=intma;
end

%Compute Exact Solution
time=0;
[qe] = exact_solution_dg(intma,coord,npoin,nelem,ngl,time,icase);
        
%Compute Initial Mass and Energy
[mass0,energy0] = compute_Mass_and_Energy(qe,intma,coord,wnq,psi,nelem,ngl,nq);

%Create Periodic BC Pointer
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

%Compute RK Time-Integration Coefficients
[a0,a1,beta] = compute_ti_coefficients(stages);

%Allocate movie arrays
dt=compute_dt(q0,intma,ngl,nelem,Courant_max,dx_min);
ntime=time_final/dt;
iplot=round(ntime/nplots);
num_plot_frames=2*nplots;
p_movie=zeros(npoin,num_plot_frames);
u_movie=zeros(npoin,num_plot_frames);
time_movie=zeros(num_plot_frames,1);
mass_movie=zeros(num_plot_frames,1);
energy_movie=zeros(num_plot_frames,1);

%Time Integration
itime=0;
while time<time_final
   %Estimate time-step
   dt=compute_dt(q0,intma,ngl,nelem,Courant_max,dx_min);
   if (dt<1e-16) 
       return
   end
   courant=Courant_max;
   time=time + dt;
   itime=itime+1;

   %Make sure we finish at Time_Final
   if time>time_final
       time=time - dt;
       dt=time_final - time;
       time=time + dt;
   end

   if mod(itime,iplot) == 0
       disp(['itime =  ',num2str(itime),' time = ', num2str(time),' dt = ', num2str(dt),' courant = ', num2str(courant)]);
   end
   
   for s=1:stages
      
      %Create RHS Matrix
      rhs = create_rhs_volume_entropy_stable(qp,intma,coord,npoin,nelem,ngl,wgl,dpsi,volume_flux_method,cgdg_method);
      rhs = create_rhs_flux_entropy_stable(rhs,qp,intma,nelem,ngl,diss,surface_flux_method,cgdg_method);
 
      %Apply Communicator
      rhs = create_rhs_communicator(rhs,space_method,intma,intma_cg,coord,Mmatrix,npoin,npoin_cg,nelem,ngl,nq,wnq,psi,0);

      %Solve System
      qp=a0(s)*q0 + a1(s)*q1 + dt*beta(s)*rhs;
                  
      %Update
      rhs_rk4(:,:,s)=rhs(:,:);
      q1=qp;
   end %ik
   
   %Do RK4 Solution
   if (stages == 4)
       qp(:,:)=q0(:,:) + dt/6.0*( rhs_rk4(:,:,1) + 2*rhs_rk4(:,:,2)...
                      + 2*rhs_rk4(:,:,3) + rhs_rk4(:,:,4) );
   end
   
   %Filter Solution
   if (mod(itime,ifilter) == 0 && ifilter > 0)
       qp = apply_filter_dg(qp,f,intma,nelem,ngl); %No DSS
       % if strcmp(space_method,'dg') %DG
       %     qp = apply_filter_dg(qp,f,intma,nelem,ngl); %No DSS
       % else
       %     qp = apply_filter_cg(qp,f,intma,coord,npoin,nelem,ngl,nq,wnq,psi,iperiodic);
       %     qp = apply_dss(qp,intma,intma_cg,coord,Mmatrix,npoin,npoin_cg,nelem,ngl,nq,wnq,psi,ivolume_integrate);
       % end
   end
   
   %Limit Solution
   limit_element=zeros(nelem,1);
   if limit == 0 || limit == 1
       [limit_element,qp,qmodal] = limiter_nodal_v3(qp,vdm,vdm_inv,intma,nelem,ngl,limit);
   elseif limit == 2
       [limit_element,qp,qmodal] = limiter_nodal_bp(qp,q0,vdm,vdm_inv,intma,nelem,ngl);
   end
   
   %DSS Solution for CG after Limiting
%    if (limit > 0 || ( mod(itime,ifilter) == 0 && ifilter > 0 ) )
%        if ( strcmp(space_method,'cgd') || strcmp(space_method,'cgc')) %CG
%           qp = apply_dss(qp,intma,intma_cg,coord,Mmatrix,npoin,npoin_cg,nelem,ngl,nq,wnq,psi,ivolume_integrate);
%        end
%    end
    % if (limit > 0)
    %    if ( strcmp(space_method,'cgd') || strcmp(space_method,'cgc')) %CG
    %       qp = apply_dss(qp,intma,intma_cg,coord,Mmatrix,npoin,npoin_cg,nelem,ngl,nq,wnq,psi,ivolume_integrate);
    %    end
    % end
   
   % %Plot Modes
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
        ylabel('u(x,t)','FontSize',18);
        set(gca, 'FontSize', 18);
        hold off;
                
%         
%         subplot(3,1,3); %Mass Conservation
%         ymax=max(mass_movie);
%         ymin=min(mass_movie);
%         tmax=max(time_movie);
%         tmin=min(time_movie);
%         tt=time_movie(1:i);
%         yt=mass_movie(1:i);
%         plot_handle=plot(tt,yt,'r-','LineWidth',2);
%         xlabel('Time','FontSize',18);
%         ylabel('\Delta M','FontSize',18);
%         set(gca, 'FontSize', 18);
% %        axis([ tmin tmax ymin ymax]);
%         hold on;
%         tt=time_movie(1:i);
%         yt=energy_movie(1:i);
%         plot_handle=plot(tt,yt,'b-','LineWidth',2);
%         legend('\Delta M','\Delta E');
%         hold on;
        
        M_i=getframe(gcf);
        M(i)=M_i;
        pause(0.01);
        hold off;
    end
    % if store_movie == 1
    %     file_movie=[main_text '_e' num2str(nelem) '_p' num2str(nop) '_q' num2str(noq) '_diss' num2str(diss) '.avi'];
    %     movie2avi(M,file_movie,'fps',5);
    % end
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
            ye(i)=qe(ip,1);
        end
        plot_handle=plot(x,y,'r-','LineWidth',2);
        hold on;
        %plot_handle=plot(x,ye,'k:','LineWidth',2);
    end
    %axis([xmin xmax hmin hmax]);
    axis([xmin xmax -2.5 2.5]);
    title_text=[main_text ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time)];
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