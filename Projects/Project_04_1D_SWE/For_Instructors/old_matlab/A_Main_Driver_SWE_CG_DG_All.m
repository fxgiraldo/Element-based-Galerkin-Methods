%---------------------------------------------------------------------%
%This code computes the 1D Shallow Water Equations using either CG or 
%DG method with either LG or LGL points for interpolation and integration.
%The time-integration is accomplished via 2nd, 3rd Order, or 3rd Order 4-stage
%RK.
%Written by F.X. Giraldo on 05/2012
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
kstages=4; %RK1, RK2 or RK3, or RK4
dt=1e-3; %time-step, fraction of one revolution
time_final=1.0; %final time in revolutions
integration_points=1; %=1 for LGL and =2 for LG
integration_type=2; %=1 is inexact and =2 is exact
space_method_type=1; %=1 for CG and =2 for DG

icase=4; %case number: 1 is a Gaussian with flat bottom; 2 is Gaussian with linear bottom
         %3 is Gaussian with Parabolic bottom; 4 is Standing Wave (linear)
         %with Analytic Solution
xmu=0.0; %filtering strength: 1 is full strength and 0 is no filter
ifilter=1; %time-step frequency that the filter is applied.
filter_type=2; %=1 is Modal Hierarchical and =2 is regular Legendre
diss=1; %=1 dissipation (Rusanov Flux) and =0 no dissipation (Central Flux)
delta_nl=0; %=0 linear and =1 nonlinear

if (icase == 4)
    delta_nl=0;
end

%Store Constants
eps=1e-15;
ntime=time_final/dt;

%Set Plotting Text for Spatial Method
if space_method_type == 1
    method_text = ['CG'];
else
    method_text = ['DG'];
end

nopp=[1 2 4 8 16 32 64];
inop_begin=1;
inop_end=5;
for inop=inop_begin:inop_end
    nop=nopp(inop);

    %Turn off filter for 1st Order Polynomials
    if nop == 1
        xmu=0;
    end

    %Compute Interpolation and Integration Points
    ngl=nop+1;
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

    icount=0;
    switch nop
        case (1)
            ibeg=32;
            iskip=32;
            iend=128;
        case (2)
            ibeg=16;
            iskip=16;
            iend=64;
        case (4)
            ibeg=8;
            iskip=8;
            iend=32;
        case (8)
            ibeg=4;
            iskip=4;
            iend=16;
        case (16)
            ibeg=2;
            iskip=2;
            iend=8;
        case (32)
            ibeg=1;
            iskip=1;
            iend=4;
         case (64)
            ibeg=1;
            iskip=1;
            iend=2;
    end %switch

    for nelem=ibeg:iskip:iend 
        disp(['nop nelem  = ',num2str(nop),'  ',num2str(nelem)]);
        icount=icount + 1;
        npoin=nop*nelem + 1;

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
        m1=0;
        for e=1:nelem
           dx=coord(ngl,e)-coord(1,e);
           jac=dx/2;
            for l=1:nq
              wq=wnq(l)*jac;
              for j=1:ngl
                 h=qe(1,j,e)+qb(j,e);
                 U=qe(2,j,e);
                 m1=m1 + wq*h*psi(j,l);
              end
            end
        end
        mass0=m1;
        
        %Create Periodic BC Pointer
        for i=1:npoin
           iperiodic(i)=i;
        end
        %iperiodic(npoin)=iperiodic(1);

        %Create Mass Matrix
        if space_method_type == 1 %CG
            Mmatrix=create_mass_cg(intma,coord,npoin,nelem,ngl,nq,wnq,psi,iperiodic);
            Mmatrix_inv=inv(Mmatrix);
        elseif space_method_type == 2 %DG
            mass=zeros(ngl,ngl,nelem);
            mass_inv=zeros(ngl,ngl,nelem);
            mass_t=zeros(ngl,ngl);
            mass = create_mass_dg(coord,nelem,ngl,nq,wnq,psi);
            for e=1:nelem
               mass_t(:,:)=mass(:,:,e);
               mass_inv(:,:,e)=inv(mass_t(:,:));
            end %ie
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

           %disp(['itime =  ',num2str(itime),' time = ', num2str(time),' courant = ', num2str(courant)]);
           for ik=1:kstages

              %Create RHS Matrix
              rhs = create_rhs(qp,qb,coord,nelem,ngl,nq,wnq,psi,dpsi,eps,gravity,delta_nl);

              %Apply Communicator
              if space_method_type == 1
                  rhs = apply_bcs_cg(rhs,qp,qb,gravity,nelem,ngl,delta_nl);
                  rhs = apply_dss(rhs,intma,Mmatrix_inv,npoin,nelem,ngl,iperiodic);
              elseif space_method_type == 2
                  rhs = create_flux(rhs,qp,qb,nelem,ngl,diss,eps,gravity,delta_nl);

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

        end %itime

        %Compute Exact Solution
        [qe,qb,gravity] = exact_solution_dg(coord,nelem,ngl,time,icase,eps);

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
              end
            end
        end
        mass_loss(icount,inop)=abs(m1-mass0);%/mass0;
            
        %Compute Norm
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
        h_norm(icount,inop)=sqrt( h_top );   
        u_norm(icount,inop)=sqrt( u_top );   
        npoin_total(icount,inop)=npoin;
    end %nelem
end %nop

h=figure;
figure(h);
for inop=inop_begin:inop_end
    switch inop
    case {1}
        plot_handle=semilogy(npoin_total(:,inop),h_norm(:,inop),'r-');
    case (2)
        plot_handle=semilogy(npoin_total(:,inop),h_norm(:,inop),'b-o');
    case (3)
        plot_handle=semilogy(npoin_total(:,inop),h_norm(:,inop),'g-x');
    case (4)
        plot_handle=semilogy(npoin_total(:,inop),h_norm(:,inop),'k-+');
    case (5)
        plot_handle=semilogy(npoin_total(:,inop),h_norm(:,inop),'m-*');
    case (6)
        plot_handle=semilogy(npoin_total(:,inop),h_norm(:,inop),'c-s');
    case (7)
        plot_handle=semilogy(npoin_total(:,inop),h_norm(:,inop),'y-d');
    end %switch
    set(plot_handle,'LineWidth',2);
    hold on
end %for
title_text=[main_text, ' , T = ' num2str(time)];
title([title_text],'FontSize',18);
xlabel('N_p','FontSize',18);
ylabel('||h||_2','FontSize',18);
legend('N=1','N=2','N=4','N=8','N=16','N=32','N=64');
set(gca, 'FontSize', 18);
axis([20 140 1e-12 1e0]);

h=figure;
figure(h);
for inop=inop_begin:inop_end
    switch inop
    case {1}
        plot_handle=semilogy(npoin_total(:,inop),u_norm(:,inop),'r-');
    case (2)
        plot_handle=semilogy(npoin_total(:,inop),u_norm(:,inop),'b-o');
    case (3)
        plot_handle=semilogy(npoin_total(:,inop),u_norm(:,inop),'g-x');
    case (4)
        plot_handle=semilogy(npoin_total(:,inop),u_norm(:,inop),'k-+');
    case (5)
        plot_handle=semilogy(npoin_total(:,inop),u_norm(:,inop),'m-*');
    case (6)
        plot_handle=semilogy(npoin_total(:,inop),u_norm(:,inop),'c-s');
    case (7)
        plot_handle=semilogy(npoin_total(:,inop),u_norm(:,inop),'y-d');
    end %switch
    set(plot_handle,'LineWidth',2);
    hold on
end %for
title_text=[main_text, ' , T = ' num2str(time)];
title([title_text],'FontSize',18);
xlabel('N_p','FontSize',18);
ylabel('||U||_2','FontSize',18);
legend('N=1','N=2','N=4','N=8','N=16','N=32','N=64');
set(gca, 'FontSize', 18);
axis([20 140 1e-12 1e0]);

h=figure;
figure(h);
for inop=inop_begin:inop_end
    switch inop
    case {1}
        plot_handle=semilogy(npoin_total(:,inop),mass_loss(:,inop),'r-');
    case (2)
        plot_handle=semilogy(npoin_total(:,inop),mass_loss(:,inop),'b-o');
    case (3)
        plot_handle=semilogy(npoin_total(:,inop),mass_loss(:,inop),'g-x');
    case (4)
        plot_handle=semilogy(npoin_total(:,inop),mass_loss(:,inop),'k-+');
    case (5)
        plot_handle=semilogy(npoin_total(:,inop),mass_loss(:,inop),'m-*');
    case (6)
        plot_handle=semilogy(npoin_total(:,inop),mass_loss(:,inop),'c-s');
    case (7)
        plot_handle=semilogy(npoin_total(:,inop),mass_loss(:,inop),'y-d');
    end %switch
    set(plot_handle,'LineWidth',2);
    hold on
end %for
title_text=[main_text, ' , T = ' num2str(time)];
title([title_text],'FontSize',18);
xlabel('N_p','FontSize',18);
ylabel('\Delta Mass','FontSize',18);
legend('N=1','N=2','N=4','N=8','N=16','N=32','N=64');
set(gca, 'FontSize', 18);
axis([20 140 1e-16 1e0]);