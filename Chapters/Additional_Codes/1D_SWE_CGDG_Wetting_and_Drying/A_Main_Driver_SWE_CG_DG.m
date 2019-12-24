%---------------------------------------------------------------------%
%This code computes the 1D Shallow Water Equations using either CG or
%DG method with either LG or LGL points for interpolation and integration.
%The time-integration is accomplished via 2nd, 3rd Order, or 3rd Order 4-stage
%RK.
%Limiting is used: Positivity Preserving (limit == 1); 
%                  Bound-Preserving (limit == 2)
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%Modified by FXG on 6/14/2014
%Code sent to Christian Sorenson for his MS Thesis 12/24/19
%---------------------------------------------------------------------%
clear all;
close all;

tic

%Input Data
nelem=100; %Number of Elements
nop=1;    %Interpolation Order

kstages=3; %RK2 or RK3
dt=0.001; %time-step for most test cases
time_final=0.5; %final time for most test cases
itime_print=100;
nplots=20; %plotting variable - Number of Frames plotted
%iplot=10;
iplot_movie=1;
store_movie=0;
integration_points=1; %=1 for LGL and =2 for LG
integration_type=2; %=1 is inexact and =2 is exact
space_method_type=2; %=1 for CG and =2 for DG

limit=1; %=0 no limiter, 
         %=1 use Shu Positivity Preserving Limiter, 
         %=2 use Bound-Preserving Limiter, 
         %=3 use both Limiters
icase=1; %case number: 
            %1 is a Gaussian with flat bottom (WORKS!!), 
            %2 is Gaussian with linear bottom (WORKS!!)
            %3 is Gaussian with Gaussian bottom (WORKS!!), 
            %4 is WD with linear bottom (works for C<= 0.01)
            %5 is Shu Test Case: WD with with still water with quadratic island (MOMENTUM
            %is TOO LARGE!!)
            %6 is WD with Parabolic Sloshing Bowl (h_eps=1e-2) (WORKS!!)
            %7 is WD with still water with linear bottom (decreasing slope:
            %WORKS!!) (h_eps=1e-3)
            %8 is WD with still water with linear bottom (increasing slope: WORKS!!)
            %9 is WD with still water with Triangular Island (WORKS!!)
            %10 is WD with still water and increasing slope quadratic beach (WORKS!!)
            %11 is WD with still water and decreasing slope quadratic beach (WORKS!!)
            %12 is still water with quadratic bathymetry
            %13 is Balzano Test 1: Linear bathymetry (seems to WORK)
            %14 is Balzano Test 2: Linear Bathymetry with Shelf (seems to WORK)
            %15 is Balzano Test 3: Linear Bathymetry with TidePool (seems to WORK)
            %16 is WD with linear Bathymetry on big domain
xmu=0.05; %filtering strength: 1 is full strength and 0 is no filter
ifilter=1; %time-step frequency that the filter is applied.
filter_type=2; %=1 is Modal Hierarchical and =2 is regular Legendre
diss=1; %=1 dissipation (Rusanov Flux) and =0 no dissipation (Central Flux)
delta_nl=1; %=0 linear and =1 nonlinear

if (icase >= 13 && icase <= 16)
    dt=1; %time-step for Balzano Tests
    time_final=10000; %final time for Balzano tests
    itime_print=1000;
end

%Store Constants
eps=1e-15;
h_eps=1e-3; %thinnest layer of water allowed 
ngl=nop + 1;
npoin=nop*nelem + 1;
ntime=time_final/dt;

%Turn off filter for 1st Order Polynomials
if nop == 1
    xmu=0;
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
[f,vdm,vdm_inv] = filter_init(ngl,xgl,xmu,filter_type);

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
mass=0;
for e=1:nelem
    %Jacobians
    dx=coord(ngl,e)-coord(1,e);
    jac=dx/2;
    for l=1:nq
        wq=wnq(l)*jac;
        for j=1:ngl
            h_k=psi(j,l);
            mass=mass + wq*(qe(1,j,e)+qb(j,e))*psi(j,l);
        end
    end
end
mass0=mass;

%Create Periodic BC Pointer
for i=1:npoin
    iperiodic(i)=i;
end
%iperiodic(npoin)=iperiodic(1);

%Create Mass Matrix
if space_method_type == 1
    [Mmatrix_h,Mmatrix_u] = create_mass_cg(intma,coord,npoin,nelem,ngl,nq,wnq,psi,iperiodic);
elseif space_method_type == 2
    mass=zeros(ngl,ngl,nelem);
    mass_inv=zeros(ngl,ngl,nelem);
    mass_t=zeros(ngl,ngl);
    rhs_p=zeros(ngl,1);
    rhs_u=zeros(ngl,1);
    rhs_inv_p=zeros(ngl,1);
    rhs_inv_u=zeros(ngl,1);
    mass = create_mass_dg(coord,nelem,ngl,nq,wnq,psi);
    for ie=1:nelem
        temp(:,:)=mass(:,:,ie);
        mass_inv(:,:,ie)=inv(temp(:,:));
    end %ie
end

j=0;
for e=1:nelem
   for i=1:ngl
       j=j+1;
       qe_p(j) = qe(1,i,e);
       qe_u(j) = qe(2,i,e);
       qe_b(j) = qb(i,e);
       x_p(j)  = coord(i,e);
   end
end
% subplot(2,1,1);
% plot(x_p,qe_p,'-sb');
% hold on
% plot(x_p,-qe_b,'-dk');
% subplot(2,1,2);
% plot(x_p,qe_u,'-sb');
% hold on
%Initialize State Vector
qp=qe; q0=qe; q1=qe;
iplot=ntime/nplots;
iframe=0;
uml=0;
counter = 0;
%Time Integration
for itime=1:ntime
    
    time=time + dt;
    u_max=-100000;
    for e=1:nelem
        for i=1:ngl
            %--------------------New precautions
            if (abs(qp(1,i,e)+qb(i,e)) < h_eps)
                u_max=max(u_max,0);
                counter=e;
            else
                u_max=max(u_max, qp(2,i,e)/(qp(1,i,e)+qb(i,e)));
                counter=e;
            end
            %--------------------New
        end
    end
    courant=u_max*dt/dx_min;
%     if (uml < u_max)
%         uml=u_max
%         counter
%     end
    
    if (mod(itime,itime_print) == 0)
        disp(['itime =  ',num2str(itime),' time = ', num2str(time),' courant = ', num2str(courant)]);
    end
    for ik=1:kstages 
        switch kstages
            case 1  %RK1
                switch ik
                    case 1
                        a0=1;
                        a1=0;
                        beta=1;
                end %ik
            case 2  %RK2
                switch ik
                    case 1
                        a0=1;
                        a1=0;
                        beta=1;
                    case (2)
                        a0=0.5;
                        a1=0.5;
                end %ik
            case 3 %RK3
                switch ik
                    case 1
                        a0=1;
                        a1=0;
                        beta=1;
                    case (2)
                        a0=3.0/4.0;
                        a1=1.0/4.0;
                        beta=1.0/4.0;
                    case (3)
                        a0=1.0/3.0;
                        a1=2.0/3.0;
                        beta=2.0/3.0;
                end %ik
        end %kstages
        dtt=dt*beta;

        %Zero Momentum for Dry Points
        %qp = zero_momentum(qp,qb,nelem,ngl,h_eps);

        %Create RHS Matrix
        rhs = create_rhs(qp,qb,coord,nelem,ngl,nq,wnq,psi,dpsi,h_eps,gravity,delta_nl);

        %Apply Communicator
        if space_method_type == 1
            rhs = apply_bcs_cg(rhs,qp,qb,gravity,nelem,ngl,delta_nl);
            rhs = apply_dss(rhs,intma,Mmatrix_h,Mmatrix_u,npoin,nelem,ngl,iperiodic);
        elseif space_method_type == 2
            rhs = create_flux(rhs,qp,qb,nelem,ngl,diss,h_eps,gravity,delta_nl,icase,time);

            %Multiply by Inverse Mass matrix
            for ie=1:nelem

                for j=1:ngl
                    for i=1:ngl
                        mass_t(i,j)=mass_inv(i,j,ie);
                    end %i
                    rhs_p(j)=rhs(1,j,ie);
                    rhs_u(j)=rhs(2,j,ie);
                end %j
                rhs(1,:,ie)=mass_t*rhs_p;
                rhs(2,:,ie)=mass_t*rhs_u;
            end
        end

        %Solve System
        qp=a0*q0 + a1*q1 + dtt*rhs;
        
         %Implement Positivity Limiter-----------------------------
         if (limit == 1 || limit == 3)
             %[qp] = limiter_Shu_TVB(qp,qb,coord,vdm,nelem,ngl);
             [qp] = limiter_Shu_Positivity_Preserving(qp,nelem,psi,ngl,nq,wnq,qb,h_eps);
             %[qp] = limiter_Shu_WD(qp,nelem,psi,ngl,nq,wnq,qb);
         end
   
        %Update
        q1=qp;
    end %ik
    
    %Filter Solution
    if (ifilter > 0)
    if (mod(itime,ifilter) == 0)
        qp = apply_filter_dg(qp,f,nelem,ngl);
        %qp=rhs;
    end
    end
    
    %Implement Bound-Preserving Limiter
    if (limit == 2 || limit == 3)
        qp = limiter_modal_bp(qp,qb,vdm,vdm_inv,nelem,ngl);
    end
       
    %Update Q
    q0=qp;
    
    if (iplot_movie == 1)
        if (mod(itime,iplot) == 0)
            iframe=iframe + 1;
            
            for e=1:nelem
                for i=1:ngl
                    p_movie(i,e,iframe)=(qp(1,i,e)+0*qb(i,e));
                    u_movie(i,e,iframe)=qp(2,i,e);
                end
            end
            time_movie(iframe)=time;
            
            %Compute Mass
            mass=0;
            for e=1:nelem
                %Jacobians
                dx=coord(ngl,e)-coord(1,e);
                jac=dx/2;
                for l=1:nq
                    wq=wnq(l)*jac;
                    for j=1:ngl
                        h_k=psi(j,l);
                        mass=mass + wq*(qp(1,j,e)+qb(j,e))*psi(j,l);
                    end
                end
            end
            mass_movie(iframe)=abs(mass-mass0)/mass0;
        end
    end

end %itime


j=0;
for e=1:nelem
   for i=1:ngl
       j=j+1;
       qp_p(j) = qp(1,i,e);
       qp_u(j) = qp(2,i,e);
       x_p(j)  = coord(i,e);
   end
end
% figure;
% subplot(2,1,1);
% plot(x_p,qp_p,'-ro');
% %title_text=[main_text ' Diss = ' num2str(diss) ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time_movie(i))];
% %title([title_text],'FontSize',18);
% 
% subplot(2,1,2);
% plot(x_p,qp_u,'-ro');


if (iplot_movie == 1)
    pmax=max(p_movie(:));
    pmin=min(-qb(:));
    %pmax=min(+1,pmax);
    xmax=max(max(coord));
    xmin=min(min(coord));
    %figure;
    figure('Position',[1 1 800 800])
    for i=1:iframe
        subplot(2,1,1);
        for e=1:nelem
            for j=1:ngl
                x(j)=coord(j,e);
                y(j)=p_movie(j,e,i);
                yb(j)=-qb(j,e);
            end
            plot_handle=plot(x,yb,'b-','LineWidth',2);
            hold on;
            plot_handle=plot(x,y,'r-','LineWidth',2);
        end
        axis([xmin xmax pmin pmax]);
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
        
        subplot(2,1,2);
        for e=1:nelem
            for j=1:ngl
                x(j)=coord(j,e);
                y(j)=u_movie(j,e,i);
            end
            plot_handle=plot(x,y,'r-');
            set(plot_handle,'LineWidth',2);
            hold on;
        end
        umax=max(max(max(u_movie)));
        umin=min(min(min(u_movie)));
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
        
%         subplot(3,1,3);
%         xt=time_movie(1:i);
%         yt=mass_movie(1:i);
%         plot_handle=plot(xt,yt,'r-');
%         set(plot_handle,'LineWidth',2);
%         
%         t1=0;
%         t2=time_movie(iframe);
%         m1=min(mass_movie);
%         m2=max(mass_movie);
%         
%         %axis([t1 t2 ]);
%         %title([title_text],'FontSize',18);
%         xlabel('Time','FontSize',18);
%         ylabel('\Delta M','FontSize',18);
%         set(gca, 'FontSize', 18);
%         hold on;
        
        M_i=getframe(gcf);
        M(i)=M_i;
        pause(0.5);
        hold off;
    end
    if store_movie == 1
        file_movie=[main_text '_e' num2str(nelem) '_p' num2str(nop) '_q' num2str(noq) '_diss' num2str(diss) '.avi'];
        v = VideoWriter(file_movie);
        v.FrameRate = 5; %set to 5 frames per second
        open(v);
        %movie2avi(M,file_movie,'fps',5);
        writeVideo(v,M);
        close(v);
    end
end

%Print Max and Mins
max(max(qp(1,:,:)))
min(min(qp(1,:,:)))
max(max(qp(2,:,:)))
min(min(qp(2,:,:)))

