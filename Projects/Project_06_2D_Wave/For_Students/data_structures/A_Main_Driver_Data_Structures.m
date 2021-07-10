%---------------------------------------------------------------------%
%A sample driver for how to use the Data Structures for Project 4
%MA4245
%Written by F.X. Giraldo on 9/17/2014
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
nelx=4; %Number of Elements
nely=2;
nop=4;    %Interpolation Order
ngl=nop+1;
noq=nop;
nq=noq+1;
plot_grid=1; %=0 don't plot, =1 plot

%Compute LGL Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);

%Compute Legendre Cardinal functions and derivatives
[psi,dpsi,xnq,wnq] = lagrange_basis(ngl,nq,xgl);

%-------------------------DATA STRUCTURES for STUDENT-------------------%
%Create Grid
[coord,intma,bsido,iperiodic,npoin,nelem,nboun,nface] = create_grid_2d(nelx,nely,nop,xgl,plot_grid);
disp(["Grid Created and stored in COORD, INTMA, and BSIDO"]);

%Compute Metric Terms
[ksi_x,ksi_y,eta_x,eta_y,jac] = metrics(coord,intma,psi,dpsi,wnq,nelem,ngl,nq);
disp(["Metric terms created and stored in KSI_X, KSI_y, ETA_x, ETA_y, and Jac"]);

%Compute Side/Edge Information
[iside,jeside] = create_side(intma,bsido,npoin,nelem,nboun,nface,ngl);
[face,imapl,imapr] = create_face(iside,intma,nface,ngl);
disp(["Face arrays created and stored in Face, IMAPL, and IMAPR"]);
[nx,ny,jac_face] = compute_normals(face,intma,coord,nface,ngl,nq,wnq,psi,dpsi);
disp(["Face normals created and stored in NX, NY, and Jac_face"]);
[face] = create_face_periodicity(iside,face,coord,nface,nboun);
%-------------------------DATA STRUCTURES for STUDENT-------------------%

%Construct Initial Conditions and do Time-integration

