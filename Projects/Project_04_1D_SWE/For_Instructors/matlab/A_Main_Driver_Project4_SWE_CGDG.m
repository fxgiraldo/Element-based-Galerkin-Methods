%---------------------------------------------------------------------%
%This code computes the 1D Shallow Water Equations using a unified CGDG 
%method using AGGP storage as described in F.X. Giraldo's text 
%"Introduction to Element-based Galerkin Methods using Tensor-Product Bases"
%
%It can use either LG or LGL points for interpolation and integration.
%The time-integration is accomplished via 2nd, 3rd, or 4th Order RK.
%
%This is the solution set to Project 4: 1D Shallow Water
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
nelem=64; %Number of Elements
nop=4;    %Interpolation Order

dt=0.002;
time_final=1; %final time in revolutions
nplots=10; %plotting variable - Number of Frames ploted
plot_movie=1;
store_movie=0;
plot_figures=0;
plot_elements=0;
plot_modes=0;
integration_points=1; %=1 for LGL and =2 for LG
integration_type=2; %=1 is inexact and =2 is exact
space_method='dg'; %cgc=CGc; cgd=CGd; dg=DG
                              
icase=4; %=1 is a Gaussian with flat bottom; 
         %2 is Gaussian with linear bottom;
         %3 is Gaussian with Parabolic bottom; 
         %4 is Standing Wave (linear) with Analytic Solution
         %5 is FXG Riemann problem;
         %6 is Simone Riemann problem
         %7 is Rupert Klein's linear multiscale problem (Ch. 7 in FXG book)
mu=1.0; %filtering strength: 1 is full strength and 0 is no filter
ifilter=0; %time-step frequency that the filter is applied.
filter_type=2; %=1 is Modal Hierarchical and =2 is regular Legendre
diss=1; %=1 dissipation (Rusanov Flux) and =0 no dissipation (Central Flux)
delta_nl=1; %=0 linear and =1 nonlinear
if (icase == 4 || icase == 7)
    delta_nl=0;
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
main_text=[space_method ':' integration_text];

%Compute Lagrange Polynomial and derivatives
[psi,dpsi] = lagrange_basis3(ngl,nq,xgl,xnq);
   
%Compute Filter Matrix
[f,vdm,vdm_inv] = filter_init(ngl,xgl,mu,filter_type);

%Create Grid
[coord,intma_cg]=create_grid(ngl,nelem,xgl,icase);
dx_min=coord(2,1)-coord(1,1);

%Form Global Matrix Pointer
npoin_cg=nop*nelem + 1;
npoin_dg=ngl*nelem;
intma=zeros(ngl,nelem);
if ( strcmp(space_method,'cgc') )
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
if (strcmp(space_method,'dg') )
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
[qe,qb,gravity] = exact_solution(intma,coord,npoin,nelem,ngl,time,icase,eps);
        
%Compute Initial Mass
m1=0; m2=0;
for e=1:nelem
   dx=coord(ngl,e)-coord(1,e);
   jac=dx/2;
    for l=1:nq
      wq=wnq(l)*jac;
      for j=1:ngl
         jp=intma(j,e);
         h=qe(jp,1)+qb(jp);
         U=qe(jp,2);
         m1=m1 + wq*h*psi(j,l);
         m2=m2 + wq*(U)*psi(j,l);
      end
    end
end
mass0=m1;
energy0=m2;

%Create Periodic BC Pointer
periodicity=zeros(npoin,1);
for i=1:npoin
   periodicity(i)=i;
end
%periodicity(npoin)=periodicity(1);

%Create Local/Element Mass
Mmatrix = create_Mmatrix(intma_cg,coord,npoin_cg,nelem,ngl,nq,wnq,psi);

%Initialize State Vector
qp=qe; q0=qe; q1=qe;
iframe=0;

%Compute RK Time-Integration Coefficients
[RKA,RKB,RKC,stages] = compute_ti_coefficients_LSRK45();
dq=zeros(npoin,2);

%Time Integration
for itime=1:ntime
   time=time + dt;
   u_max=-100000;
   for e=1:nelem
       for i=1:ngl
           ip=intma(i,e);
           u_max=max(u_max, qp(ip,2)/(qp(ip,1)+qb(ip)));
       end
   end
   courant=u_max*dt/dx_min;
    
   if (mod(itime,100) == 0)
       disp(['itime =  ',num2str(itime),' time = ', num2str(time),' courant = ', num2str(courant)]);
   end
   
   %LSRK45 Stages
   for s = 1:stages
       %Create RHS vector: volume, flux, and communicator
       R = create_rhs_volume(qp,qb,intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,gravity,delta_nl);
       R = create_rhs_flux(R,qp,qb,intma,nelem,ngl,diss,gravity,delta_nl);
       R = create_rhs_communicator(R,space_method,intma,intma_cg,coord,Mmatrix,npoin,npoin_cg,nelem,ngl,nq,wnq,psi,periodicity,0);
       %Solve System
       for I=1:npoin
           dq(I,:) = RKA(s)*dq(I,:) + dt*R(I,:);
           qp(I,:) = qp(I,:) + RKB(s)*dq(I,:);
       end
       if periodicity(npoin) == periodicity(1)
           qp(npoin,:)=qp(1,:); %periodicity
       end
   end %s
    
   %Filter Solution
   if (mod(itime,ifilter) == 0 && ifilter > 0)
      qp = apply_filter(qp,f,intma,nelem,ngl);
   end
   
   %DSS Solution for CG
   if (mod(itime,ifilter) == 0 && ifilter > 0 )
       if ( strcmp(space_method,'cgd') || strcmp(space_method,'cgc')) %CG
          qp = apply_dss(qp,intma,intma_cg,coord,Mmatrix,npoin,npoin_cg,nelem,ngl,nq,wnq,psi,periodicity,1);
       end
   end
   
   %Plot Modes
   if (plot_modes == 1 && mod(itime,iplot) == 0)
       [test] = plot_modes(qmodal,Ipointer,limit_element,intma,coord,qp,nelem,ngl,nq,nop,noq,time,diss,main_text);           
   end

   %Update Q
   q0=qp;
   
   if (plot_movie == 1 && mod(itime,iplot) == 0)
      [p_movie,u_movie,time_movie,mass_movie,energy_movie,L2_norm,iframe] = store_movie_frames(qp,q0,qb,intma,coord,wnq,psi,npoin,nelem,ngl,nq,p_movie,u_movie,time_movie,mass_movie,energy_movie,iframe,mass0,energy0,time,icase,eps);           
   end 
   
end %itime

if (plot_movie == 1)
    
    figure;
    for i=1:iframe
        subplot(2,2,1); %H Solution
        hmax=max(p_movie(:));
        hmin=min(-qb);
        xmax=max(coord(:));
        xmin=min(coord(:));
        for e=1:nelem
            for j=1:ngl
                jp=intma(j,e);
                x(j)=coord(j,e);
                y(j)=p_movie(jp,i);
                yb(j)=-qb(jp);
            end
            plot_handle=plot(x,y,'r-','LineWidth',2);
            hold on;
            plot_handle=plot(x,yb,'b-','LineWidth',2);
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
        
        subplot(2,2,2); %U Solution
        umax=max(u_movie(:));
        umin=min(u_movie(:));
        for e=1:nelem
            for j=1:ngl
                jp=intma(j,e);
                x(j)=coord(j,e);
                y(j)=u_movie(jp,i);
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
            i1=intma(1,e);
            i2=intma(ngl,e);
            x1=coord(1,e);
            x2=coord(ngl,e);
            y1=u_movie(i1,i);
            y2=u_movie(i2,i);
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
end

if (plot_figures)
    
    %Exact Solution
    [qe,qb,gravity] = exact_solution(intma,coord,npoin,nelem,ngl,time,icase,eps);
    
    figure; %H Plot
    hmax=max(qp(:,1));
    hmin=min(qp(:,1));
    xmax=max(max(coord(:,:)));
    xmin=min(min(coord(:,:)));
    for e=1:nelem
        for i=1:ngl
            ip=intma(i,e);
            x(i)=coord(i,e);
            y(i)=qp(ip,1);
            ye(i)=qe(ip,1);
            yb(i)=-qb(ip);
        end
        plot(x,y,'r-','LineWidth',2);
        hold on;
        plot(x,ye,'k:','LineWidth',2);
        plot(x,yb,'b-','LineWidth',2);
    end
    axis([xmin xmax hmin hmax]);
    title_text=[main_text ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time)];
    title([title_text],'FontSize',18);
    hold on;

    %Plot Elements
    if (plot_elements == 1)
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
    ylabel('h','FontSize',18);
    set(gca, 'FontSize', 18);
    hold off;
    
    figure; %U Plot
    hmax=max(qp(:,2));
    hmin=min(qp(:,2));
    xmax=max(max(coord(:,:)));
    xmin=min(min(coord(:,:)));
    for e=1:nelem
        for i=1:ngl
            ip=intma(i,e);
            x(i)=coord(i,e);
            y(i)=qp(ip,2);
            ye(i)=qe(ip,2);
        end
        plot(x,y,'r-','LineWidth',2);
        hold on;
        plot(x,ye,'k:','LineWidth',2);
    end
    axis([xmin xmax hmin hmax]);
    title_text=[main_text ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time)];
    title([title_text],'FontSize',18);
    hold on;
    
    %Plot Elements
    if (plot_elements == 1)
        for e=1:nelem
            i1=intma(1,e);
            i2=intma(ngl,e);
            x1=coord(1,e);
            x2=coord(ngl,e);
            y1=qp(i1,2);
            y2=qp(i2,2);
            plot(x1,y1,'ko');
            plot(x2,y2,'ko');
        end
    end

    title([title_text],'FontSize',18);
    xlabel('x','FontSize',18);
    ylabel('U','FontSize',18);
    set(gca, 'FontSize', 18);
    hold off;
end