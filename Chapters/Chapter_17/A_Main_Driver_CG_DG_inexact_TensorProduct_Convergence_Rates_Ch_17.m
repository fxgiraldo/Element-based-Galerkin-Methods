%---------------------------------------------------------------------%
%This code computes the 2D Advection Equation using the CG/DG methods
%with 2nd Order RK and tensor product of 1D basis function with 
%Inexact Integration (Inexact Integration and Tensor-Product basis
%functions)
%Written by F.X. Giraldo on 6/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

space_method='cg'; %CG or DG
kstages=3;  %2=RK2, 3=RK3
dt=1; %time-step, Changes automatically to keep Courant_max fixed!
Courant_max=0.25;
time_final=1; %final time in revolutions
nplots=10; %Number of Frames in movie
store_movie=0;
plot_grid=0; %=0 don't plot, =1 plot %DO NOT CHANGE!

icase=1; %case number: 1 is a Gaussian in CW, 2 is Gaussian along X, 
         %3 is Gaussian along Y, 4 is Gaussian along Diagonal, 
         %5 is square wave in CW, 6 is square wave along X.
xmu=0.05; %filtering strength: 1 is full strength and 0 is no filter
ifilter=0; %time-step frequency that the filter is applied.

%Store Constants
if (icase ==1)
    c=2*pi;
elseif (icase ==2)
    c=1;
elseif (icase ==3)
    c=1;
elseif (icase ==4)
    c=1;
elseif (icase ==5)
    c=2*pi;
elseif (icase ==6)
    c=1;
end
%ntime=time_final/dt;
dt=dt*c;
time_final=time_final*c;

nopp=[1 2 4 8 16];
inop_begin=1;
inop_end=5;
for inop=inop_begin:inop_end
nop=nopp(inop);
ngl=nop + 1;
nq=ngl;

icount=0;
switch nop
    case (1)
        ibeg=8;
        iskip=8;
        iend=48;
    case (2)
        ibeg=4;
        iskip=4;
        iend=24;
    case (4)
        ibeg=2;
        iskip=2;
        iend=12;
    case (8)
        ibeg=1;
        iskip=1;
        iend=6;
    case (16)
        ibeg=1;
        iskip=1;
        iend=3;
end %switch
main_text=[space_method ':'];

%Compute LGL Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);

%Compute Legendre Cardinal functions and derivatives
[psi,dpsi,xnq,wnq] = lagrange_basis(ngl,nq,xgl);

%Compute Filter Matrix
f = filter_init(ngl,xgl,xmu);

for nel=ibeg:iskip:iend
icount=icount + 1;
t0=cputime;

nelx=nel;
nely=nel;
ntime=time_final/dt;

%-------------------------DATA STRUCTURES for STUDENT-------------------%
%Create Grid
[coord,intma,bsido,iperiodic,npoin,nelem,nboun,nface] = create_grid_2d(nelx,nely,nop,xgl,plot_grid);
                            
%Compute Metric Terms
[ksi_x,ksi_y,eta_x,eta_y,jac] = metrics(coord,intma,psi,dpsi,wnq,nelem,ngl,nq);

%Compute Side/Edge Information
[iside,jeside] = create_side(intma,bsido,npoin,nelem,nboun,nface,ngl);
[face,imapl,imapr] = create_face(iside,intma,nface,ngl);
[nx,ny,jac_face] = compute_normals(face,intma,coord,...
                   nface,ngl,nq,wnq,psi,dpsi);
face=create_face_periodicity(iside,face,coord,nface,nboun);
%-------------------------DATA STRUCTURES for STUDENT-------------------%
     
%Compute Exact Solution
time=0;

[qa,ua,va]=exact_solution(coord,npoin,time,icase);

q0=zeros(nelem,ngl,ngl);
qe=zeros(nelem,ngl,ngl);
ue=zeros(nelem,ngl,ngl);
ve=zeros(nelem,ngl,ngl);

% Reshape Arrays
for e=1:nelem
    for j=1:ngl
    for i=1:ngl
        ip=intma(e,i,j);
        qe(e,i,j)=qa(ip);
        ue(e,i,j)=ua(ip);
        ve(e,i,j)=va(ip);
    end %i
    end %j
end %e

%Create Mass Matrix
Mmatrix_inv=zeros(nelem,ngl,ngl);
Mmatrix = create_Mmatrix2d_TensorProduct_inexact(jac,intma,nelem,ngl);  
if (space_method == 'cg')
    Mmatrix = apply_dss(Mmatrix,intma,iperiodic,ngl,npoin,nelem);
end
%Invert Mass
for e=1:nelem
    Mmatrix_inv(e,:,:)=1/Mmatrix(e,:,:);
end

%Compute Courant Number
[Courant,vel,ds,dt] = compute_Courant(ua,va,intma,coord,nelem,ngl,dt,Courant_max);
ntime=round(time_final/dt);
dt=time_final/ntime;
Courant=vel*dt/ds;
%pause(3);

%Initialize State Vector
q1=qe;
q0=qe;
qp=qe;
iplot=round(ntime/nplots);
iframe=0;

%Time Integration
for itime=1:ntime
   itime;
   time=time + dt;
   timec=time/(c);
   timec;

    for ik=1:kstages
        switch kstages
            case 2  %RK2
                switch ik
                    case 1
                        a0=1;
                        a1=0;
                        beta=1;
                    case (2)
                        a0=0.5;
                        a1=0.5;
                        beta=0.5;
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
      
      %Construct RHS
      rhs = compute_rhs_TensorProduct_inexact(qp,ue,ve,ksi_x,ksi_y,eta_x,eta_y,jac,...
	        dpsi,nelem,ngl);
        
      %Construct Communicator: DSS for CG or Fluxes for DG
      if (space_method == 'cg')
          rhs = apply_dss(rhs,intma,iperiodic,ngl,npoin,nelem);
      elseif (space_method == 'dg')
        rhs = compute_flux_TensorProduct_inexact(rhs,qp,ue,ve,face,nx,ny,jac_face,...
            nface,ngl,imapl,imapr);      
      end %if
      rhs=rhs./Mmatrix; 
      
      %Evolve forward in Time
      for e=1:nelem
          qp(e,:,:)=a0*q0(e,:,:) + a1*q1(e,:,:) + dtt*rhs(e,:,:);
      end

      %Filter Solution
       if (ifilter > 0)
       if (mod(itime,ifilter) == 0)
          rhs = apply_filter2D_dg(qp,f,nelem,ngl);
          if (space_method == 'cg')
               rhs=rhs.*jac;
               rhs = apply_dss(rhs,intma,iperiodic,ngl,npoin,nelem);
               rhs=rhs./Mmatrix; 
          end
          qp=rhs;
       end
       end

      %Update
      q1=qp;
      
   end %ik

   
   %Update Q
   q0=qp;
     
end %itime

%Compute Exact Solution
[qa,ua,va] = exact_solution(coord,npoin,time,icase);

% Reshape Arrays
for e=1:nelem
    for j=1:ngl
    for i=1:ngl
        ip=intma(e,i,j);
        qe(e,i,j)=qa(ip);
        ue(e,i,j)=ua(ip);
        ve(e,i,j)=va(ip);
    end %i
    end %j
end %e

%Compute Norm
top=0;
bot=0;
for e=1:nelem
    for j=1:ngl
    for i=1:ngl
        top=top + (q0(e,i,j)-qe(e,i,j))^2;
        bot=bot + qe(e,i,j)^2;
    end %i
    end %j
end %e
l2_norm=sqrt( top/bot );
l2_norm_total(icount,inop)=l2_norm;
npoin_total(icount,inop)=npoin;
if space_method == 'dg'
    npoin_total(icount,inop)=nelem*ngl^2;
end
t1=cputime;
wallclock_time(icount,inop)=t1-t0;
disp([' nop  = ' num2str(nop),' nel = ' num2str(nel),' l2 = ',num2str(l2_norm),' cpu = ' num2str(t1-t0) ]);

end %nel
end %nopp

h=figure;
figure(h);
for inop=inop_begin:inop_end
    switch inop
    case (1)
        plot_handle=semilogy(npoin_total(:,inop),l2_norm_total(:,inop),'r-');
    case (2)
        plot_handle=semilogy(npoin_total(:,inop),l2_norm_total(:,inop),'b-o');
    case (3)
        plot_handle=semilogy(npoin_total(:,inop),l2_norm_total(:,inop),'g-x');
    case (4)
        plot_handle=semilogy(npoin_total(:,inop),l2_norm_total(:,inop),'k-+');
    case (5)
        plot_handle=semilogy(npoin_total(:,inop),l2_norm_total(:,inop),'m-*');
    case (6)
        plot_handle=semilogy(npoin_total(:,inop),l2_norm_total(:,inop),'c-s');
end %switch
set(plot_handle,'LineWidth',2);
hold on
end
title_text=[main_text ' RK', num2str(kstages), ' T = ' num2str(time)];
title([title_text],'FontSize',18);
xlabel('N_P','FontSize',18);
ylabel('Normalized L^2 Error','FontSize',18);
legend('N=1','N=2','N=4','N=8','N=16');
set(gca, 'FontSize', 18);
axis([0 10000 1e-4 1e0]);

h=figure;
figure(h);
for inop=inop_begin:inop_end
    switch inop
    case (1)
        plot_handle=semilogy(wallclock_time(:,inop),l2_norm_total(:,inop),'r-');
    case (2)
        plot_handle=semilogy(wallclock_time(:,inop),l2_norm_total(:,inop),'b-o');
    case (3)
        plot_handle=semilogy(wallclock_time(:,inop),l2_norm_total(:,inop),'g-x');
    case (4)
        plot_handle=semilogy(wallclock_time(:,inop),l2_norm_total(:,inop),'k-+');
    case (5)
        plot_handle=semilogy(wallclock_time(:,inop),l2_norm_total(:,inop),'m-*');
    case (6)
        plot_handle=semilogy(wallclock_time(:,inop),l2_norm_total(:,inop),'c-s');
end %switch
set(plot_handle,'LineWidth',2);
hold on
end
title_text=[main_text ' RK', num2str(kstages), ' T = ' num2str(time)];
title([title_text],'FontSize',18);
xlabel('Wallclock Time (S)','FontSize',18);
ylabel('Normalized L^2 Error','FontSize',18);
legend('N=1','N=2','N=4','N=8','N=16');
set(gca, 'FontSize', 18);

