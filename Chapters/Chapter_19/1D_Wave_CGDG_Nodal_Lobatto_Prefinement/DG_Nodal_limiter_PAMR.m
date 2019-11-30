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
nelem=30; %Number of Elements
nop=4;    %Interpolation Order
noq=nop; %NOP and NOQ need to be the same
nop_min=1; %minimum interpolation order
nop_max=4; %maximum interpolation order
p_tol=1e-4; %relative change to trigger refinement (norm2(high)/norm2(total))
p_tol2=1e-4; %max norm of element solution
adapt=1; %1=true, 0=false
if (nop ~= nop)
    disp(['Error: NOQ must be equal to NOP']);
    noq=nop;
end

kstages=3; %RK2 or RK3
dt=0.01; %time-step, fraction of one revolution
Courant_max=0.25;
time_final=0.25; %final time in revolutions
nplots=10; %plotting variable - Number of Frames ploted
iplot_movie=1;
iplot_solution=1; %Switch to Plot or Not.
iplot_modes=1;
store_movie=0;
store_plot=0;

limit=1; %=0 no limiting, =1 yes to limiting
% limiter_type=3; %=1 Krividonoa; =2 FVM; 3=Bound-Preserving Limiter; 4=TVB Limiter
% m_limiter=1e+2;
icase=2; %case number: 1 is a Gaussian, 2 is a square wave, 3 is
         %Gaussian with source and 4 is square wave with source, 5 is a
         %step function and 6 is a gaussian and a square wave
xmu=0.0; %filtering strength: 1 is full strength and 0 is no filter
ifilter=1; %time-step frequency that the filter is applied.
diss=1; %=1 dissipation (Rusanov Flux) and =0 no dissipation (Central Flux)

%Modify according to ADAPT flag
if (adapt == 0)
    nop_min=nop;
    nop_max=nop;
end

%Store Constants
ngl=nop + 1;
nq=noq + 1;
npoin=nop*nelem + 1;
ntime=time_final/dt;

element_order=zeros(nelem,1);
element_order(1:nelem)=nop;
%Compute i,e -> I pointer
I=0;
for e=1:nelem
    for i=1:ngl
        I=I+1;
        Ipointer(i,e)=I;
    end
end

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
dt=Courant_max*ds/u;
ntime=round(time_final/dt)
dt=time_final/ntime
Courant=u*dt/ds

%Compute Exact Solution
time=0;
qe = exact_solution_dg(coord,nelem,ngl,time,icase);
fe = source_function_dg(coord,nelem,ngl,time,icase);

%Create Mass Matrix
mass=zeros(ngl,1);
jac = create_jacobian_dg(coord,nelem,ngl);

%Initialize State Vector
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
   
   disp(['itime =  ',num2str(itime),' time = ', num2str(time),' courant = ', num2str(Courant)]);
   for ik=1:kstages      
      %Create RHS Matrix
      rhs = create_rhs_dg(qp,jac,nelem,ngl,wgl_matrix,psi_matrix,dpsi_matrix,u,diss,element_order);

      %Solve
      for e=1:nelem
          m=element_order(e);
          n=m+1;
          mass(1:n)=wgl_matrix(1:n,m)*jac(e);
          qp(1:n,e)=alpha(ik,1)*q0(1:n,e) + alpha(ik,2)*q1(1:n,e) + dt*beta(ik)*rhs(1:n,e)./mass(1:n);
      end %ie
      
      %Update
      q1=qp;
   end %ik
   
   %Filter Solution
   if (mod(itime,ifilter) == 0)
      rhs = apply_filter_dg(qp,f,nelem,ngl);
      qp=rhs;
   end
   
   %P-Adaptivity
   if adapt == 1
      [qp,element_order] = change_polynomial_order(qp,VDM_Matrix,VDM_Inv_Matrix,nelem,element_order,p_tol,p_tol2,nop_min,nop_max);
   end %if

   %Limit Solution
    if limit == 1
        qp = limiter_Shu_Positivity_Preserving(qp,nelem,element_order,nop_max,wgl_matrix,Interpolation_Matrix);
    end
      
   %Update Q
   q0=qp;
   
   %Map Solution to Quadrature Points
   qq=zeros(ngl,nelem);
   for e=1:nelem
        n=element_order(e)+1;
        for k=1:ngl
            for i=1:n
                qq(k,e)=qq(k,e) + Interpolation_Matrix(k,i,n-1)*qp(i,e);
            end
        end
   end

   %Map to Modal Solution
   qmodal=zeros(ngl,nelem);
   for e=1:nelem
       N=element_order(e)+1;
       i=N;   
       %Get Modal Coefficients
       qmodal(1:i,e)=VDM_Inv_Matrix(1:i,1:i,i-1)*qp(1:i,e);
   end
    
   if (iplot_modes == 1)
       
       II=Ipointer';
       QQ=qmodal';
       icheck=[1:ngl];
       
       %-------PLOT MODES------%
       subplot(3,1,1);
       %bar(II(1:nelem,icheck),QQ(1:nelem,icheck),'stacked');
       bar(QQ(1:nelem,icheck));
       %xlabel('Methods','FontSize',18);
       %set(gca,'XTickLabel',{'1-norm','2-norm','inf-norm','fro-norm'});    
       %ylabel('Mode Values','FontSize',18);
       title_text=[ 'Modes ' ];
       %title([title_text],'FontSize',18);
       %legend('0','1','2','3','4','5','6','7','8');
       title_text=['DG Nodal: Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time)];
       title([title_text],'FontSize',18);
       qmin=min(min(qmodal));
       qmax=max(max(qmodal));
       axis([ 1 nelem qmin qmax]);
       set(gca,'FontSize',18);
       %hold on;
       
       %Mark Off Element Boundaries
%        for e=1:nelem
%         x1=Ipointer(1,e);
%         x2=Ipointer(ngl,e);
%         y1=0;
%         y2=0;
%         plot(x1,y1,'k.','LineWidth',2);
%         plot(x2,y2,'k.','LineWidth',2);
%        end

%        %-------PLOT Limiter------%
%        subplot(3,1,2);
%        bar(limit_element);
%        qmin=min(limit_element);
%        qmax=max(limit_element);
%        qmin=0;
%        qmax=1;
%        axis([ 1 nelem qmin qmax]);
%        set(gca,'FontSize',18);
       
       %-------PLOT Element Order------%
       subplot(3,1,2);
       bar(element_order);
       qmin=min(element_order);
       qmax=max(element_order);        
       qmin=0;
       qmax=nop_max;
       axis([ 1 nelem qmin qmax]);
       set(gca,'FontSize',18);

       %-------PLOT Solution------%
       subplot(3,1,3)
       for e=1:nelem
            for i=1:ngl
                x(i)=coord(i,e);
                y(i)=qq(i,e);
            end
            plot_handle=plot(x,y,'r-');
            set(plot_handle,'LineWidth',2);
            hold on;
       end %e
        axis([ -1 +1 -0.25 1.25]);
        set(gca,'FontSize',18);
%        M_i=getframe(gcf);
%        M_modes(i)=M_i;
        pause(0.1);
       hold off;
       
   end %IPLOT_MODES
   
   if (iplot_movie == 1)
      if (mod(itime,iplot) == 0)
           iframe=iframe + 1;
           qi_movie(:,:,iframe)=qq;
           time_movie(iframe)=time;
      end
   end 
   
end %itime

if (iplot_movie == 1)
    figure;
    for i=1:iframe
        
        for e=1:nelem
            for j=1:ngl
                x(j)=coord(j,e);
                y(j)=qi_movie(j,e,i);
            end
            plot_handle=plot(x,y,'r-');
            set(plot_handle,'LineWidth',2);
            hold on;
        end
        axis([-1 +1 -0.25 +1.25]);
        title_text=['DG LGL: Diss = ' num2str(diss) ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time_movie(i))];
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
        file_movie=['DG_LGL_n' num2str(nelem) '_p' num2str(nop) '_q' num2str(noq) '_diss' num2str(diss) '.avi'];
        movie2avi(M,file_movie,'fps',5);
    end
end

%Compute Exact Solution
qe = exact_solution_dg(coord,nelem,ngl,time,icase);

%Compute Norm
top=0;
bot=0;

for ie=1:nelem
   for i=1:ngl
	   top=top + (qq(i,ie)-qe(i,ie))^2;
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
    ip=0;
    for ie=1:nelem
    for i=1:ngl
          ip=ip+1;
          q_sol(ip)=q_sol(ip) + qq(i,ie);
          qe_sol(ip)=qe_sol(ip) + qe(i,ie);
          x_sol(ip)=coord(i,ie);
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
    for ie=1:nelem
        x1=coord(1,ie);
        x2=coord(ngl,ie);
        y1=qq(1,ie);
        y2=qq(ngl,ie);
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

    if (diss == 0)
       file_ps=['dg_lgl_n' num2str(nelem) 'p' num2str(nop)];
        legend('DG LGL','Exact');	
    elseif (diss == 1)
       file_ps=['dg_lgl_upwind_n' num2str(nelem) 'p' num2str(nop)];
       legend('DG LGL','Exact');	
    end

    title_text=['DG Nodal: Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm) ', T = ' num2str(time)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
    
    if store_plot == 1
        eval(['print ' file_ps ' -depsc']);
    end
end


dt
ds
Courant
q_max=max(max(qq))
q_min=min(min(qq))
l2_norm

for e=1:nelem
    i=element_order(e);
    ii=i+1;
    for k=1:ii
        if (qp(k,e) < 0)
            disp(['e = ',num2str(e),' k = ',num2str(k),' qp = ',num2str(qp(k,e))]);
        end
    end
end
