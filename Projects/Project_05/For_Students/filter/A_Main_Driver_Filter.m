%---------------------------------------------------------------------%
%A sample driver for how to use the Filter routines for Project 4
%MA4245
%Written by F.X. Giraldo on 9/17/2014
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

%Input Data
ngl=4;    %Interpolation Order
space_method='cg'; %CG or DG
xmu=0.05; %filtering strength: 1 is full strength and 0 is no filter
ifilter=1; %time-step frequency that the filter is applied. ifilter=1 means filter every time-step

dt=2*pi/100;
time_final=2*pi;
ntime=time_final/dt;

%Compute LGL Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);

%Compute Filter Matrix
f = filter_init(ngl,xgl,xmu);

%BUILD GRID, Initial Conditions....

%Time Integration
for itime=1:ntime
   time=time + dt;
       
  %Evolve forward in Time
  for e=1:nelem
      qp(e,:,:)=a0*q0(e,:,:) + a1*q1(e,:,:) + dtt*rhs(e,:,:);
  end

  %Apply Filter to Solution
   if (mod(itime,ifilter) == 0)
      rhs = apply_filter2D_dg(qp,f,nelem,ngl);
      if (space_method == 'cg') %for CG we need to DSS the solution
           rhs=rhs.*jac;
           rhs = apply_dss(rhs,intma,iperiodic,ngl,npoin,nelem);
           rhs=rhs./Mmatrix; 
      end
      qp=rhs; %for DG nothing needs to be done since discontinuous solution is OK
   end
end %itime

