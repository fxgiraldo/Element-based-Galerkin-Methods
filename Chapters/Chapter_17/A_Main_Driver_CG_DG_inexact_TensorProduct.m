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

tic

%Input Data
nel=4; %Number of Elements
nop=4;    %Interpolation Order
noq=nop; %DO NOT CHANGE!
space_method='dg'; %CG or DG
kstages=3;  %2=RK2, 3=RK3
dt=1; %time-step, Changes automatically to keep Courant_max fixed!
Courant_max=0.25;
time_final=0.25; %final time in revolutions
nframes=10; %Number of Frames in movie
store_movie=0;
plot_grid=1; %=0 don't plot, =1 plot
grid_rotation_angle=45;

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
nelx=nel;
nely=nel;
ngl=nop + 1;
nq=noq + 1;
ntime=time_final/dt;

%Compute LGL Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);

%Compute Legendre Cardinal functions and derivatives
[psi,dpsi,xnq,wnq] = lagrange_basis(ngl,nq,xgl);

%Compute Filter Matrix
f = filter_init(ngl,xgl,xmu);

%-------------------------DATA STRUCTURES for STUDENT-------------------%
%Create Grid
[coord,intma,bsido,iperiodic,npoin,nelem,nboun,nface] = create_grid_2d(nelx,nely,nop,xgl,plot_grid);

%Rotate Grid
[coord_rotated] = rotate_grid_v2(coord,intma,npoin,nelem,ngl,plot_grid,grid_rotation_angle);

%Compute Metric Terms
[ksi_x,ksi_y,eta_x,eta_y,jac] = metrics(coord_rotated,intma,psi,dpsi,wnq,nelem,ngl,nq);

%Compute Side/Edge Information
if (space_method == 'dg')
    [iside,jeside] = create_side(intma,bsido,npoin,nelem,nboun,nface,ngl);
    [face,imapl,imapr] = create_face(iside,intma,nface,ngl);
    [nx,ny,jac_face] = compute_normals(face,intma,coord,nface,ngl,nq,wnq,psi,dpsi);
    [face] = create_face_periodicity(iside,face,coord,nface,nboun);
    [nx_rotated,ny_rotated] = rotate_normals(nx,ny,nface,nq,grid_rotation_angle);
    nx=nx_rotated;
    ny=ny_rotated;
end

%Store Rotated COORDS
coord=coord_rotated;
%-------------------------DATA STRUCTURES for STUDENT-------------------%

%Compute Exact Solution
time=0;
[qa,ua,va]=exact_solution(coord,npoin,time,icase);
q0=zeros(nelem,ngl,ngl);
qe=zeros(nelem,ngl,ngl);
ue=zeros(nelem,ngl,ngl);
ve=zeros(nelem,ngl,ngl);

%Reshape Arrays
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
Mmatrix = create_Mmatrix2d_TensorProduct_inexact(jac,intma,nelem,ngl);  
if (space_method == 'cg')
    Mmatrix = apply_dss(Mmatrix,intma,iperiodic,ngl,npoin,nelem);
end
%Invert Mass
for e=1:nelem
    Mmatrix_inv(e,:,:)=1/Mmatrix(e,:,:);
end

%Compute Courant Number
[Courant,vel,ds,dt] = compute_Courant(ua,va,intma,coord,nelem,ngl,dt,Courant_max)
ntime=round(time_final/dt);
dt=time_final/ntime;
Courant=vel*dt/ds;
%pause(3);

%Initialize State Vector
q1=qe;
q0=qe;
qp=qe;
% iplot=round(ntime/nframes);
% iframe=0;

%Intialize Movie
[qi_movie,time_movie,xi,yi,iframe] = initialize_movie(coord,nframes);

%Time Integration
for itime=1:ntime
   itime;
   time=time + dt;
   timec=time/(c);
   timec;
   
   disp(['itime time courant = ',num2str(itime),' ',num2str(timec),' ',num2str(Courant)]);

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
       if (mod(itime,ifilter) == 0)
          rhs = apply_filter2D_dg(qp,f,nelem,ngl);
          if (space_method == 'cg')
               rhs=rhs.*jac;
               rhs = apply_dss(rhs,intma,iperiodic,ngl,npoin,nelem);
               rhs=rhs./Mmatrix; 
          end
          qp=rhs;
       end

      %Update
      q1=qp;
      
   end %ik

   
   %Update Q
   q0=qp;
   
   %PLOT Solution for MOVIES
   [qi_movie,time_movie,iframe] = update_movie(qi_movie,time_movie,q0,coord,intma,npoin,nelem,ngl,ntime,itime,nframes,iframe,xi,yi,timec);

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

%Compute gridpoint solution
q_sol=zeros(npoin,1);
qe_sol=zeros(npoin,1);
lhowm=zeros(npoin,1);
for ie=1:nelem
    for j=1:ngl
    for i=1:ngl
      ip=intma(ie,i,j);
      lhowm(ip)=lhowm(ip)+1;
      q_sol(ip)=q_sol(ip) + q0(ie,i,j);
      qe_sol(ip)=qe_sol(ip) + qe(ie,i,j);
    end %i
    end %j
end
for i=1:npoin
   q_sol(i)=q_sol(i)/lhowm(i);
   qe_sol(i)=qe_sol(i)/lhowm(i);
end

%Plot Solution
h=figure;
figure(h);
qi=griddata(coord(:,1),coord(:,2),q_sol,xi,yi,'cubic');
[cl,h]=contourf(xi,yi,qi);
colorbar('SouthOutside');
xlabel('X','FontSize',18);
ylabel('Y','FontSize',18);
axis image
title_text=[space_method ' : Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm)];
title([title_text],'FontSize',18);
set(gca, 'FontSize', 18);
%file_ps=['se_n' num2str(nelem) 'p' num2str(nop)];
%eval(['print ' file_ps ' -depsc']);

%Plot Movie
Movie_frames = plot_movie(qi_movie,time_movie,xi,yi,iframe);
if (store_movie == 1)
    movie2avi(M,'dg_inexact_TensorProduct_movie.avi','fps',5);
end

nop
nelem
dt=dt/(c);
dt
Courant
q_max=max(max(max(q0)))
q_min=min(min(min(q0)))
l2_norm

toc
