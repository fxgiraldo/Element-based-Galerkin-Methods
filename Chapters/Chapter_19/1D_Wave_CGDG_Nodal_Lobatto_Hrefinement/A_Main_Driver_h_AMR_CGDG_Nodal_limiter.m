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
space_method='dg'; %CG or DG
nelem=10; %Number of Elements
nop=4;    %Interpolation Order
noq=nop; %NOP and NOQ need to be the same
nop_min=1; %minimum interpolation order
nop_max=4; %maximum interpolation order
refine_tol=1e-2; %relative change to trigger refinement (norm2(high)/norm2(total))
coarsen_tol=1e-4; %max norm of element solution
padapt=0; %1=true, 0=false
hadapt=1; %1=true, 0=false
ref_level_max=2; %number of refinement levels
h_eps=0.34; %for ad-hoc initial adaptation

if (nop ~= nop)
    disp(['Error: NOQ must be equal to NOP']);
    noq=nop;
end

kstages=3; %RK2 or RK3
dt=0.001; %time-step, fraction of one revolution
Courant_max=0.25;
time_final=1.0; %final time in revolutions
nplots=10; %plotting variable - Number of Frames ploted
iplot_movie=0;
iplot_solution=0; %Switch to Plot or Not.
iplot_modes=1;
store_movie=0;
store_plot=0;

diss=1; %=1 dissipation (Rusanov Flux) and =0 no dissipation (Central Flux)
limit=0; %=0 no limiting, =1 yes to limiting
% limiter_type=3; %=1 Krividonova; =2 FVM; 3=Bound-Preserving Limiter; 4=TVB Limiter
% m_limiter=1e+2;
icase=1; %case number: 1 is a Gaussian, 2 is a square wave, 3 is
         %Gaussian with source and 4 is square wave with source, 5 is a
         %step function and 6 is a gaussian and a square wave
         
%Filtering Not Yet implemented for CG
xmu=0.05; %filtering strength: 1 is full strength and 0 is no filter
ifilter=1; %time-step frequency that the filter is applied.

%Modify according to ADAPT flag
if (padapt == 0)
    nop_min=nop;
    nop_max=nop;
end

%Store Constants
ngl=nop + 1;
nq=noq + 1;
npoin=nop*nelem + 1;
ntime=time_final/dt;

%Compute LGL Points
xgl_matrix=zeros(ngl,nop);
wgl_matrix=zeros(ngl,nop);
for i=1:nop
    n=i+1;
    xgl=zeros(n,1);
    wgl=zeros(n,1);
    [xgl,wgl]=legendre_gauss_lobatto(n);    
    xgl_matrix(1:n,i)=xgl(1:n);
    wgl_matrix(1:n,i)=wgl(1:n);
end
xgl=xgl_matrix(1:ngl,nop);
wgl=wgl_matrix(1:ngl,nop);
xnq=xgl;
wnq=wgl;

%Construct Interpolation Matrices: VDM
VDM_Matrix=zeros(ngl,ngl,nop);
VDM_Inv_Matrix=zeros(ngl,ngl,nop);
for i=1:nop
    ii=i+1;
    xi=zeros(ii,1);
    xi(1:ii)=xgl_matrix(1:ii,i);
    L = legendre_basis_modal_new(ii,ii,xi);
    VDM_Matrix(1:ii,1:ii,i)=L(1:ii,1:ii);
    temp=inv(L);
    VDM_Inv_Matrix(1:ii,1:ii,i)=temp(1:ii,1:ii);
end

%Construct Projection Matrix from Low to High
Interpolation_Matrix=zeros(ngl,ngl,nop);
for i=1:nop
    ii=i+1;
    xi=zeros(ii,1);
    xi(1:ii)=xgl_matrix(1:ii,i);
    xk(1:ngl)=xgl_matrix(1:ngl,nop);
    [psi,dpsi] = lagrange_basis3(ii,ngl,xi,xk);
    for m=1:ii
        for n=1:ngl
            Interpolation_Matrix(n,m,i)=psi(m,n); %maps I -> K
        end
    end
end

%Compute Lagrange Polynomial and derivatives
psi_matrix=zeros(ngl,ngl,nop);
dpsi_matrix=zeros(ngl,ngl,nop);
for i=1:nop
    n=i+1;
    xe=zeros(n,1);
    psi=zeros(n,n);
    dpsi=zeros(n,n);
    xe(1:n)=xgl_matrix(1:n,i);
    [psi,dpsi] = lagrange_basis3(n,n,xe,xe);
    psi_matrix(1:n,1:n,i)=psi(1:n,1:n);
    dpsi_matrix(1:n,1:n,i)=dpsi(1:n,1:n);
end
psi(1:n,1:n)=psi_matrix(1:n,1:n,nop);
dpsi(1:n,1:n)=dpsi_matrix(1:n,1:n,nop);

%Compute Filter Matrix
[f,vdm,vdm_inv] = filter_init(ngl,xgl,xmu);

%Create Grid
[coord,intma]=create_grid_dg(ngl,nelem,xgl);
ds=coord(2,1)-coord(1,1);
u=2;

%H-adapt Arrays
[parent,children,active,ref_level,elem_lev_pointer,coord,intma,nelem,nelem0]=create_hadapt_arrays(intma,coord,nelem,ngl,xgl,ref_level_max);
[sfc,nsfc]=create_space_filling_curve(nelem,nelem0,children,active);
[P1g,P2g,P1s,P2s] = L2_Projection_1D_fxg(ngl,nq);

%Apply Ad Hoc H-AMR
% if hadapt == 1
%     active = apply_hadapt_init(children,active,elem_lev_pointer,coord,ngl,hadapt_lev,h_eps);
%     [sfc,nsfc]=create_space_filling_curve(nelem,nelem0,children,active);
% end

%Compute Courant Number
[Courant,ds] = compute_Courant_number(coord,ngl,u,sfc,nsfc,dt);
% dt=Courant_max*ds/u
ntime=round(time_final/dt)
% dt=time_final/ntime
Courant=u*dt/ds;

%P-adapt Array
element_order=zeros(nelem,1);
element_order(1:nelem)=nop;

%Compute Exact Solution
time=0;
qe = exact_solution_dg(coord,nelem,ngl,time,icase);
fe = source_function_dg(coord,nelem,ngl,time,icase);

%Create Mass Matrix
jac = create_jacobian_dg(coord,nelem,ngl);
mass = create_mass_element(jac,wgl_matrix,nelem,ngl,element_order);
if strcmp(space_method,'cg')
   mass = apply_dss(mass,intma,npoin,nelem,element_order,sfc,nsfc);
end

qp=qe; q0=qe; q1=qe;
iplot=round(ntime/nplots);
%iplot=1;
iframe=0;

%Form Butcher Tableau for RK3
alpha=zeros(3,3);
beta=zeros(3,1);
alpha(1,1)=1; alpha(1,2)=0;
alpha(2,1)=3.0/4.0; alpha(2,2)=1.0/4.0;
alpha(3,1)=1.0/3.0; alpha(3,2)=2.0/3.0;
beta(1)=1; beta(2)=1.0/4.0; beta(3)=2.0/3.0;

%Time Integration
for itime=1:ntime
   time=time + dt;
   
   disp(['itime=',num2str(itime),' time=', num2str(time),' Courant=', num2str(Courant),' nelem=',num2str(nelem),' nsfc=',num2str(nsfc)]);
   for ik=1:kstages      
      %Create RHS Matrix
      rhs = create_rhs_dg(qp,jac,nelem,ngl,wgl_matrix,psi_matrix,dpsi_matrix,u,diss,element_order,space_method,sfc,nsfc);
      
      %DSS for CG
      if strcmp(space_method,'cg')
          rhs = apply_dss(rhs,intma,npoin,nelem,element_order,sfc,nsfc);
      end
      
      %Solve
      for ee=1:nsfc
          e=sfc(ee);         
          i=element_order(e);
          ii=i+1;
          qp(1:ii,e)=alpha(ik,1)*q0(1:ii,e) + alpha(ik,2)*q1(1:ii,e) + dt*beta(ik)*rhs(1:ii,e)./mass(1:ii,e);
      end %e
      
      %Update
      q1=qp;
   end %ik
   
   %Filter Solution
   if (mod(itime,ifilter) == 0)
      rhs = apply_filter_dg(qp,f,nelem,ngl,sfc,nsfc);
%       if strcmp(space_method,'cg')
%           rhs = apply_dss(rhs,intma,npoin,nelem,element_order);
%           for e=1:nelem
%               i=element_order(e);
%               ii=i+1;
%               rhs(1:ii,e)=rhs(1:ii,e)./mass(1:ii,e);
%           end %e
%       end
      qp = rhs;
   end
   
%    %P-Adaptivity
%    if padapt == 1
%       [qp,element_order] = change_polynomial_order(qp,VDM_Matrix,VDM_Inv_Matrix,nelem,element_order,p_tol,p_tol2,nop_min,nop_max);
%       mass = create_mass_element(jac,wgl_matrix,nelem,ngl,element_order);
%       if strcmp(space_method,'cg')
%           mass = apply_dss(mass,intma,npoin,nelem,element_order,sfc,nsfc);
%       end
%    end
   %H-Adaptivity
    if hadapt == 1
        [qp,active] = apply_hadapt_dynamic(qp,active,ref_level,children,parent,ref_level_max,sfc,nsfc,VDM_Inv_Matrix,nelem,ngl,refine_tol,coarsen_tol,P1g,P2g,P1s,P2s);
        [sfc,nsfc]=create_space_filling_curve(nelem,nelem0,children,active);
        [Courant,ds] = compute_Courant_number(coord,ngl,u,sfc,nsfc,dt);
    end

   %Limit Solution
    if limit == 1
        qp = limiter_Shu_Positivity_Preserving(qp,nelem,element_order,nop_max,wgl_matrix,Interpolation_Matrix);
    end
      
   %Update Q
   q0=qp;
   
   %Map to Modal Solution
   qmodal=zeros(ngl,nsfc);
   ref_level_frame=zeros(nsfc,1);
   midpoint_frame=zeros(nsfc,1);
   for ee=1:nsfc
       e=sfc(ee);
       N=element_order(e)+1;
       i=N;   
       %Get Modal Coefficients
       qmodal(1:i,ee)=VDM_Inv_Matrix(1:i,1:i,i-1)*qp(1:i,e);
       ref_level_frame(ee)=ref_level(e);
       midpoint_frame(ee)=0.5*( coord(N,e)+coord(1,e) );
   end
    
   if (iplot_modes == 1)
       
       QQ=qmodal';
       icheck=[1:ngl];
       
       %-------PLOT MODES------%
       subplot(3,1,1);
       bar(QQ(1:nsfc,icheck));
       title_text=[ 'Modes' ];
       title_text=[space_method ': Ne = ' num2str(nsfc) ', N = ' num2str(nop) ', Time = ' num2str(time)];
       title([title_text],'FontSize',18);
       qmin=min(qmodal(:));
       qmax=max(qmodal(:));
       axis([ 1 nsfc qmin qmax]);
       %xlabel('x','FontSize',18);
       ylabel('Modes','FontSize',18);
       set(gca,'FontSize',18);
             
       %-------PLOT Refinement Level------%
       subplot(3,1,2);
%        bar(midpoint_frame,ref_level_frame);
       bar(ref_level_frame);
       qmax=ref_level_max;
       qmin=0;
%        axis([ -1 +1 qmin qmax]);
       axis([ 1 nsfc qmin qmax]);
       set(gca,'FontSize',18);
       %xlabel('x','FontSize',18);
       ylabel('H-Level','FontSize',18);
       
       %-------PLOT Solution------%
       subplot(3,1,3)
       for ee=1:nsfc
            e=sfc(ee);
            for i=1:ngl
                x(i)=coord(i,e);
                y(i)=qp(i,e);
            end
            plot_handle=plot(x,y,'r-');
            set(plot_handle,'LineWidth',2);
            hold on;
       end %e
       axis([ -1 +1 -0.25 1.25]);
       set(gca,'FontSize',18);
       xlabel('x','FontSize',18);
       ylabel('q(x,t)','FontSize',18);
       
        %Plot Discontinuous Elements
        for ee=1:nsfc
            e=sfc(ee);
            x1=coord(1,e);
            x2=coord(ngl,e);
            y1=qp(1,e);
            y2=qp(ngl,e);
            plot(x1,y1,'ko');
            plot(x2,y2,'ko');
        end

%        M_i=getframe(gcf);
%        M_modes(i)=M_i;
       pause(0.1);
       hold off;
       
   end %IPLOT_MODES
   
   if (iplot_movie == 1)
      if (mod(itime,iplot) == 0)
           iframe=iframe + 1;
           qi_movie(:,:,iframe)=qp;
           time_movie(iframe)=time;
           sfc_movie(:,iframe)=sfc(:);
           nsfc_movie(iframe)=nsfc;
      end
   end 
   
end %itime

if (iplot_movie == 1)
    figure;
    for i=1:iframe      
        for ee=1:nsfc_movie(i)
            e=sfc_movie(ee,i);
            for j=1:ngl
                x(j)=coord(j,e);
                y(j)=qi_movie(j,e,i);
            end
            plot_handle=plot(x,y,'r-');
            set(plot_handle,'LineWidth',2);
            hold on;
        end
        axis([-1 +1 -0.25 +1.25]);
        title_text=[space_method ': Ne = ' num2str(nsfc_movie(i)) ', N = ' num2str(nop) ', Time = ' num2str(time_movie(i))];
        title([title_text],'FontSize',18);
        hold on;
        
        %Plot Elements
        for ee=1:nsfc_movie(i)
            e=sfc_movie(ee,i);
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
        file_movie=[space_method ':' num2str(nsfc) '_p' num2str(nop) '_q' num2str(noq) '_diss' num2str(diss) '.avi'];
        movie2avi(M,file_movie,'fps',5);
    end
end

%Compute Exact Solution
qe = exact_solution_dg(coord,nelem,ngl,time,icase);

%Compute Norm
top=0;
bot=0;

for ee=1:nsfc
   e=sfc(ee);
   for i=1:ngl
       top=top + (qp(i,e)-qe(i,e))^2;
       bot=bot + qe(i,e)^2;
   end %i
end %ie
l2_norm=sqrt( top/bot );

%Plot Solution
if (iplot_solution == 1)
    %Compute a gridpoint solution
    np=nsfc*ngl;
    q_sol=zeros(np,1);
    qe_sol=zeros(np,1);
    x_sol=zeros(np,1);
    ip=0;
    for ee=1:nsfc
        e=sfc(ee);
        for i=1:ngl
              ip=ip+1;
              q_sol(ip)=q_sol(ip) + qp(i,e);
              qe_sol(ip)=qe_sol(ip) + qe(i,e);
              x_sol(ip)=coord(i,e);
        end 
    end
    
    h=figure;
    figure(h);
    plot_handle=plot(x_sol,q_sol,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(x_sol,qe_sol,'b--');
    set(plot_handle,'LineWidth',2);
    axis([-1 +1 -0.25 +1.25]);

    %Plot Discontinuous Elements
    for ee=1:nsfc
        e=sfc(ee);
        x1=coord(1,e);
        x2=coord(ngl,e);
        y1=qp(1,e);
        y2=qp(ngl,e);
        plot(x1,y1,'ko');
        plot(x2,y2,'ko');
    end

%     %Plot Continuous Elements
%     for ie=1:nelem
%         i1=intma(1,ie);
%         i2=intma(ngl,ie);
%         x1=x_sol(i1);
%         x2=x_sol(i2);
%         y1=q_sol(i1);
%         y2=q_sol(i2);
%         plot(x1,y1,'ko');
%         plot(x2,y2,'ko');
%     end

    xlabel('x','FontSize',18);
    ylabel('q(x,t)','FontSize',18);

    if (diss == 0)
       file_ps=[space_method ':' num2str(nsfc) 'p' num2str(nop)];
        legend('DG LGL','Exact');	
    elseif (diss == 1)
       file_ps=['dg_lgl_upwind_n' num2str(nsfc) 'p' num2str(nop)];
       legend('DG LGL','Exact');	
    end

    title_text=[space_method ': Ne = ' num2str(nsfc) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm) ', T = ' num2str(time)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
    
    if store_plot == 1
        eval(['print ' file_ps ' -depsc']);
    end
end


dt
ds
Courant
q_max=max(max(qp))
q_min=min(min(qp))
l2_norm

for ee=1:nsfc
    e=sfc(ee);
    i=element_order(e);
    ii=i+1;
    for k=1:ii
        if (qp(k,e) < 0)
            disp(['e = ',num2str(e),' k = ',num2str(k),' qp = ',num2str(qp(k,e))]);
        end
    end
end
