%---------------------------------------------------------------------%
%This code computes the 1D Advection Equation using CG
%with LGL/LGR points for interpolation and integration.
%The time-integration is accomplished via 2nd, 3rd Order, or 3rd Order 4-stage
%RK.
%Written by F.X. Giraldo on 11/2023
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

addpath('fxg') %New LGL and LGR points and bases

%------------------------------------Input Data
nelem_LGL=20; %Number of Elements
nelem_LGR=1;
nop_LGL=4;    %Polynomial Order
nop_LGR=20;
nelem=nelem_LGL+nelem_LGR;
lgr_scale=1; %1=scale and 0=do not scale
lgr_scale_factor=1; %factor to scale from xgr_scale
basis_method_LGL_Only=0; %----------------------------Make sure it is 0!
form_method='strong';
xmin=-1;
xmax=+1;

kstages=4; %RK1, RK2, RK3, or RK4
dt=0.0005;
time_final=1; %final time in revolutions
nplots=100; %plotting variable - How many DT to plot
iplot_movie=1;
store_movie=0;
iplot_figures=1;
iplot_elements=0;
iplot_u_max=1;
icase=1; %=1 is a Gaussian with flat bottom; 
%------------------------------------Input Data

%Store Constants
eps=1e-15;
ntime=time_final/dt;
iplot=nplots;

%Store vectors
time_vector=zeros(ntime,1);
umax_vector=zeros(ntime,1);

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

%Compute Interpolation and Integration Points
[xgl,wgl]=legendre_gauss_lobatto(ngl_LGL);
if basis_method_LGL_Only == 1
    [xgr,wgr]=legendre_gauss_lobatto(ngl_LGR);
else
    %[xgr,wgr]=laguerre_gauss_radau_eigenvalue(ngl_LGR,true); %true=scale weights by exp(x)
    [xgr,wgr]=laguerre_gauss_radau_eigenvalue_fxg(ngl_LGR);
end

%Scale LGR element to be closer to LGL elements
dx_lgl=(xmax-xmin)/nelem_LGL;
dx_lgr=xgr(ngl_LGR)-xgr(1);
xgr_scale=1.0;
if (lgr_scale==1)
    xgr_scale=dx_lgl/dx_lgr*lgr_scale_factor;
end

%Compute Lagrange Polynomial and derivatives
[psil,dpsil] = lagrange_basis3(ngl_LGL,ngl_LGL,xgl,xgl);

if basis_method_LGL_Only == 1
    [psir,dpsir] = lagrange_basis3(ngl_LGR,ngl_LGR,xgr,xgr);
else
    [psir,dpsir,~,~] = lagrange_basis_laguerre_v3(ngl_LGR,ngl_LGR,xgr); 
end
dpsi=zeros(ngl(nelem),ngl(nelem),nelem);
for e=1:nelem_LGL
    for j=1:ngl(e)
        for i=1:ngl(e)
            dpsi(i,j,e)=dpsil(i,j);
        end
    end
end
for e=nelem_LGL+1:nelem
    for j=1:ngl(e)
        for i=1:ngl(e)
            dpsi(i,j,e)=dpsir(i,j);
        end
    end
end

%Create Grid
[coord,coord_cg,intma,jac,wnq,npoin]=create_grid(ngl,nelem,nelem_LGL,xgl,xgr,wgl,wgr,icase,xgr_scale,basis_method_LGL_Only,xmin,xmax);
dx_min=min(coord(2,1)-coord(1,1),coord(2,nelem)-coord(1,nelem));

%print LGL and LGR element sizes
e=1;
dx_LGL=coord(ngl(e),e) - coord(1,e);
e=nelem;
dx_LGR=coord(ngl(e),e) - coord(1,e);
disp([' dx_LGL = ', num2str(dx_LGL),' dx_LGR = ', num2str(dx_LGR)]);
pause(1.0);

% for e=1:nelem
%     dx=coord(ngl(e),e) - coord(1,e);
%     disp(['e =  ',num2str(e),' dx = ', num2str(dx)]);
% end

%Compute Exact Solution
time=0;
h_movie=zeros(npoin,nplots);
u_movie=zeros(npoin,nplots);
time_movie=zeros(nplots,1);
mass_movie=zeros(nplots,1);
energy_movie=zeros(nplots,1);
[qe,u0] = exact_solution(intma,coord,npoin,nelem,nelem_LGL,ngl);
        
%Compute Initial Mass and Energy
[mass0] = compute_Mass_and_Energy(qe,intma,jac,wnq,nelem,ngl);

%Create Local/Element Mass
Mmatrix = create_Mmatrix(intma,jac,wnq,npoin,nelem,ngl);

%Damping Sponge Layer (see Benacchio and Bonaventura, 2013)
%Klemp-Lily sponge
zt = max(coord_cg);
zd = coord(ngl_LGL,nelem_LGL);
z = max(coord_cg-zd,0);
gammamax = 1;
gamma = gammamax.*sin(pi/2 .* z./(zt -zd)).^2;
Igamma=1.0 - gamma;
disp(['ztop =  ',num2str(zt),' zd = ', num2str(zd)]);

%Initialize State Vector
qp=qe; q0=qe; q1=qe;
iframe=0;

%Initialize Arrays
rhs=zeros(npoin,1);
rhs_rk4=zeros(npoin,4);

%Compute RK Time-Integration Coefficients
[a0,a1,beta] = compute_ti_coefficients(kstages);

%Time Integration
time_vector(1)=time;
u_max_vector(1)=max(qp(:));
for itime=1:ntime
   time=time + dt;
   u_max=-100000;
   I_max=0;
   for e=1:nelem
       for i=1:ngl(e)
           I=intma(i,e);
           if (abs(qp(I)) > u_max)
               u_max=abs(qp(I));
               I_max=I;
           end
       end
   end
   courant=u_max*dt/dx_min;
   time_vector(itime)=time;
   u_max_vector(itime)=u_max;

   if mod(itime,iplot) == 0
       disp(['itime =  ',num2str(itime),' time = ', num2str(time),' u_max = ', num2str(u_max) ' I_max = ', num2str(I_max)]);
   end
   
   for ik=1:kstages
      %Create RHS Matrix
      rhs = create_rhs_volume(qp,u0,intma,jac,npoin,nelem,ngl,wnq,dpsi,form_method);
      
      %Apply Communicator
      rhs = create_rhs_communicator(rhs,Mmatrix,npoin);
      
      %Solve System
      for I=1:npoin
            qp(I)=( a0(ik)*q0(I) + a1(ik)*q1(I) + dt*beta(ik)*rhs(I) )*Igamma(I);
      end

      %BC
      qp(1)=0; qp(npoin)=0;

      %Update
      rhs_rk4(:,ik)=rhs(:);
      q1=qp;
   end %ik
   
   %Do RK4 Solution
   if (kstages == 4)
       for I=1:npoin
           qp(I)=(q0(I) + dt/6.0*( rhs_rk4(I,1) + 2*rhs_rk4(I,2)...
                       + 2*rhs_rk4(I,3) + rhs_rk4(I,4) ))*Igamma(I);
       end
       %BC
       qp(1)=0; qp(npoin)=0;
   end
   
   %Update Q
   q0=qp;
   if (iplot_movie == 1 && mod(itime,iplot) == 0)
      iframe=iframe+1;
      h_movie(:,iframe)=qp(:);
      time_movie(iframe)=time;
   end 
   
end %itime

if (iplot_movie == 1)
    
    hmax=max(h_movie(:));
    hmin=min(h_movie(:));
    xmax=max(coord_cg(:));
    xmin=min(coord_cg(:));
    figure;
    for i=1:iframe
        plot_handle=plot(coord_cg(:),h_movie(:,i),'r-','LineWidth',2);
        axis([xmin xmax hmin hmax]);
        title_text=['Ne = ' num2str(nelem) ', N-LGL = ' num2str(nop_LGL) ', N-LGR = ' num2str(nop_LGR) ', Time = ' num2str(time_movie(i))];
        title([title_text],'FontSize',18);
        hold on;
         
        %Plot Elements
        if (iplot_elements == 1)
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
        plot(xx,yy,'b:','LineWidth',2.0);
        
        title([title_text],'FontSize',18);
        xlabel('x','FontSize',18);
        ylabel('q(x,t)','FontSize',18);
        set(gca, 'FontSize', 18);
        hold off;
                       
        M_i=getframe(gcf);
        M(i)=M_i;
        pause(0.01);
        hold off;
    end
    % if store_movie == 1
    %     file_movie=['Advection_e=' num2str(nelem) '_p=' num2str(nop) '.avi'];
    %     movie2avi(M,file_movie,'fps',5);
    % end
end

if (iplot_figures)
        
    figure; %H Plot
    hmax=max(qp(:));
    hmin=min(qp(:));
    xmax=max(coord(:));
    xmin=min(coord(:));
    for e=1:nelem
        for i=1:ngl(e)
            I=intma(i,e);
            x(i)=coord(i,e);
            y(i)=qp(I);
            ye(i)=qe(I);
        end
        plot_handle=plot(x,y,'r-','LineWidth',2);
        hold on;
        %plot_handle=plot(x,ye,'k:','LineWidth',2);
    end
    axis([xmin xmax hmin hmax]);
    title_text=[' Ne = ' num2str(nelem) ', N-LGL = ' num2str(nop_LGL) ', N-LGR = ' num2str(nop_LGR) ', Time = ' num2str(time)];
    title([title_text],'FontSize',18);
    hold on;

    %Plot Elements
    if (iplot_elements == 1)
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
    plot(xx,yy,'b:','LineWidth',2.0);

    title([title_text],'FontSize',18);
    xlabel('x','FontSize',18);
    ylabel('q(x,t)','FontSize',18);
    set(gca, 'FontSize', 18);
    hold off;
            
end

if (iplot_u_max==1)
    figure;
    hmax=max(u_max_vector(:));
    hmin=min(u_max_vector(:));
    xmax=max(time_vector(:));
    xmin=min(time_vector(:));
    plot_handle=semilogy(time_vector,u_max_vector,'r-','LineWidth',2);
    axis([xmin xmax hmin hmax]);
    title_text=[' Ne = ' num2str(nelem) ', N-LGL = ' num2str(nop_LGL) ', N-LGR = ' num2str(nop_LGR)];
    title([title_text],'FontSize',18);
    hold on;

    xlabel('time','FontSize',18);
    ylabel('max(|q(x,t^n)|)','FontSize',18);
    set(gca, 'FontSize', 18);
end