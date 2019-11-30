%---------------------------------------------------------------------%
%This code computes the 1D Advection Equation using the DG method with
%Modal functions
%with 2nd or 3rd Order RK with Legendre Points.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
nelem=10; %Number of Elements
nop=8;    %Interpolation Order
nop_min=1; %minimum interpolation order
nop_max=8; %maximum interpolation order
p_tol_c=1e-8; %relative change to trigger refinement (norm2(high)/norm2(total))
p_tol_r=1e-4; %relative change to trigger refinement (norm2(high)/norm2(total))
p_tol2=1e-4; %max norm of element solution
adapt=1; %1=true, 0=false

kstages=3; %RK2 or RK3
dt=0.01; %time-step, fraction of one revolution
Courant_max=0.25;
time_final=1.0; %final time in revolutions
nplots=10; %plotting variable - Number of Frames ploted
iplot_movie=0;
iplot_solution=0; %Switch to Plot or Not.
iplot_filter=0;
iplot_modes=1;
store_movie=0;
store_plot=0;

limit=0; %=0 no limiting, =1 yes to limiting
icase=1; %case number: 1 is a Gaussian, 2 is a square wave, 3 is
         %Gaussian with source and 4 is square wave with source, 5 is a
         %step function and 6 is a gaussian and a square wave
xmu=0; %filtering strength: 1 is full strength and 0 is no filter
ifilter=1; %time-step frequency that the filter is applied.
diss=1; %=1 for dissipation or =0 for no dissipation

%Modify according to ADAPT flag
if (adapt == 0)
    nop_min=nop;
    nop_max=nop;
end

%Store Constants
ngl=nop + 1;
noq=nop_max;
nq=noq+1;
ntime=time_final/dt;
npoin=nq*nelem;
iplot=round(ntime/nplots);

%Compute i,e -> I pointer and Element_order array
I=0;
element_order=zeros(nelem,1);
element_order(1:nelem)=nop;
for e=1:nelem
    for i=1:ngl
        I=I+1;
        Ipointer(i,e)=I;
    end
end

%Compute LGL Points
%[xnq,wnq]=legendre_gauss(nq);
[xnq,wnq]=legendre_gauss_lobatto(nq);

%Construct Modal Basis Functions
[L,dL] = legendre_basis_modal(nq,nq,xnq); %Form the Legendre Basis up to the maximum order

%Sample Lagrange Polynomial at Element Interface (Xi=+/-1)
ns=2;
xs(1)=-1; xs(2)=+1;
% [Ls,dLs] = legendre_basis_modal(ngl,ns,xs);
[Ls,dLs] = legendre_basis_modal(nq,ns,xs);

%Compute Filter Matrix
f = filter_init_modal(ngl,xmu,iplot_filter);
% pause

%Create Grid
[coord,intma,dx]=create_grid_dg_lg(nq,nelem,xnq);
ds=coord(2,1)-coord(1,1);
%ds=coord(nq,1)-coord(1,1);
u=2;
dt=Courant_max*ds/u;
ntime=round(time_final/dt)
dt=time_final/ntime
Courant=u*dt/ds

%Compute Exact Solution in Physical Space
time=0;
qe = exact_solution_dg(coord,nelem,nq,time,icase);
fe = source_function_dg(coord,nelem,nq,time,icase);

%Get Initial Condition in Modal Space
q0=zeros(nq,nelem);
for e=1:nelem
    for i=1:ngl
        sum=0;
        for k=1:nq 
            sum=sum + wnq(k)*L(i,k)*qe(k,e);
        end
        a=2/(2*(i-1)+1);
        q0(i,e)=sum/a;
    end
end
    
%Create Mass Matrix
mass = create_mass_dg_modal(nelem,nq,nq,wnq,L,dx);
mass_inv=zeros(nq,nq,nelem);
for e=1:nelem
   temp(:,:)=mass(:,:,e);
   mass_inv(:,:,e)=inv(temp(:,:));
end %e

%Initialize State Vector
qq=zeros(nq,1);
qp=q0; q1=q0;
iplot=round(ntime/nplots);
iframe=0;

%Initialize Arrays
rhs_t=zeros(ngl,1);

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
      rhs = create_rhs_dg_modal(qp,nelem,nq,wnq,L,Ls,dL,u,diss,dx,element_order);
      
      %Solve
      qp(:)=0;
      for e=1:nelem
         N=element_order(e)+1;
         qp(1:N,e)=alpha(ik,1)*q0(1:N,e) + alpha(ik,2)*q1(1:N,e) + dt*beta(ik)*mass_inv(1:N,1:N,e)*rhs(1:N,e);
      end %e
      
      %Update
      q1=qp;
   end %ik
   
   %Filter Solution
%    if (mod(itime,ifilter) == 0)
%        for e=1:nelem
%            qp(:,e)= f(:,:)*qp(:,e);
%        end
%    end
   
   
    %P-Adaptivity
    if adapt == 1
        [qp,element_order] = change_polynomial_order_v2(qp,nelem,element_order,p_tol_c,p_tol_r,p_tol2,nop_min,nop_max);
    end %if
    
    %Limit Solution
    if limit == 1
        qp = limiter_Shu_Positivity_Preserving_v2(qp,nelem,L,nq,wnq,element_order);
    end

   %Update Q
   q0=qp;
   
   if (iplot_modes == 1)
       II=Ipointer';
       QQ=qp';
       icheck=[1:ngl];
       
       %-------PLOT MODES------%
       subplot(3,1,1);
       %bar(II(1:nelem,icheck),QQ(1:nelem,icheck),'stacked');
       bar(QQ(1:nelem,icheck));
       %xlabel('Methods','FontSize',18);
       %set(gca,'XTickLabel',{'1-norm','2-norm','inf-norm','fro-norm'});    
       ylabel('Modes','FontSize',18);
       title_text=[ 'Modes ' ];
       %title([title_text],'FontSize',18);
       %legend('0','1','2','3','4','5','6','7','8');
       title_text=['DG Modal: Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time)];
       title([title_text],'FontSize',18);
       qmin=min(min(qp));
       qmax=max(max(qp));
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

       %-------PLOT Limiter------%
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
       ylabel('P-Level','FontSize',18);
       qmin=min(element_order);
       qmax=max(element_order);        
       qmin=0;
       qmax=nop_max;
       axis([ 1 nelem qmin qmax]);
       set(gca,'FontSize',18);

       %-------PLOT Solution------%
       subplot(3,1,3)
       qq=zeros(nq,nelem);
        for ie=1:nelem
            for i=1:nq
                for k=1:ngl
                    qq(i,ie)=qq(i,ie) + L(k,i)*qp(k,ie);
                end
                x(i)=coord(i,ie);
                y(i)=qq(i,ie);
            end
            
            plot_handle=plot(x,y,'r-');
            set(plot_handle,'LineWidth',2);
            hold on;
        end
        axis([ -1 +1 -0.25 1.25]);
        set(gca,'FontSize',18);
%        M_i=getframe(gcf);
%        M_modes(i)=M_i;
        xlabel('x','FontSize',18);
        ylabel('q(x,t)','FontSize',18);
       
        %Plot Discontinuous Elements
        for e=1:nelem
            x1=coord(1,e);
            x2=coord(ngl,e);
            y1=qq(1,e);
            y2=qq(ngl,e);
            plot(x1,y1,'ko');
            plot(x2,y2,'ko');
        end
        
        qmin=min(element_order);
       
        pause(0.1);
        hold off;
       
   end
   
   if (iplot_movie == 1)
      if (mod(itime,iplot) == 0)
          iframe=iframe + 1;
          qq=zeros(nq,nelem);
          for e=1:nelem
              N=element_order(e)+1;
              for i=1:nq
                for k=1:N
                    qq(i,e)=qq(i,e) + L(k,i)*qp(k,e);
                end
                x(i)=coord(i,e);
                y(i)=qq(i,e);
            end
          end
          plot_handle=plot(x,y,'r-');
          set(plot_handle,'LineWidth',2);
          hold on;
          qi_movie(:,:,iframe)=qq;
          time_movie(iframe)=time;
      end
   end 
   
end %itime

if (iplot_movie == 1)
    figure;
    for i=1:iframe
        
        for e=1:nelem
            for j=1:nq
                x(j)=coord(j,e);
                y(j)=qi_movie(j,e,i);
            end
            plot_handle=plot(x,y,'r-');
            set(plot_handle,'LineWidth',2);
            hold on;
        end
        axis([-1 +1 -0.25 +1.25]);
        title_text=['DG Modal: Diss = ' num2str(diss) ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time_movie(i))];
        title([title_text],'FontSize',18);
        hold on;
        
        %Plot Elements
        for e=1:nelem
            x1=coord(1,e);
            x2=coord(nq,e);
            y1=qi_movie(1,e,i);
            y2=qi_movie(nq,e,i);
            plot(x1,y1,'k.','LineWidth',2);
            plot(x2,y2,'k.','LineWidth',2);
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
        file_movie=['DG_Modal_n' num2str(nelem) '_p' num2str(nop) '_q' num2str(noq) '_diss' num2str(diss) '.avi'];
        movie2avi(M,file_movie,'fps',5);
    end
end

%Compute Exact Solution
qe = exact_solution_dg(coord,nelem,nq,time,icase);

%Map Solution to Physical Space
qq=zeros(nq,nelem);
for e=1:nelem
    N=element_order(e)+1;
    for i=1:nq
        for k=1:N
            qq(i,e)=qq(i,e) + L(k,i)*q0(k,e);
        end
    end
end

%Compute Norm
top=0;
bot=0;
for e=1:nelem 
   for i=1:nq
	   top=top + (qq(i,e)-qe(i,e))^2;
       bot=bot + qe(i,e)^2;
   end %i
end %ie
l2_norm=sqrt( top/bot );

%Plot Solution
if (iplot_solution == 1)
    %Compute a gridpoint solution
    %Plot Initial Condition
    np=nelem*nq;
    q_sol=zeros(np,1);
    qe_sol=zeros(np,1);
    x_sol=zeros(np,1);
    ip=0;
    for e=1:nelem
        for i=1:nq
              ip=ip+1;
              q_sol(ip)=q_sol(ip) + qq(i,e);
              qe_sol(ip)=qe_sol(ip) + qe(i,e);
              x_sol(ip)=coord(i,e);
        end 
    end
    figure;
    plot_handle=plot(x_sol,q_sol,'r-');
    set(plot_handle,'LineWidth',2);
    
    hold on
    plot_handle=plot(x_sol,qe_sol,'b--');
    set(plot_handle,'LineWidth',2);
    axis([-1 +1 -0.25 +1.25]);

    %Plot Discontinuous Elements
    for e=1:nelem
        x1=coord(1,e);
        x2=coord(nq,e);
        y1=qq(1,e);
        y2=qq(nq,e);
        plot(x1,y1,'ko');
        plot(x2,y2,'ko');
    end

    xlabel('x','FontSize',18);
    ylabel('q(x,t)','FontSize',18);

    if (diss == 0)
       file_ps=['dg_modal_n' num2str(nelem) 'p' num2str(nop)];
        legend('DG Modal','Exact');	
    elseif (diss == 1)
       file_ps=['dg_modal_upwind_n' num2str(nelem) 'p' num2str(nop)];
       legend('DG Modal','Exact');	
    end

    title_text=['DG Modal: Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm) ', T = ' num2str(time)];
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
nop_min
nop_max


