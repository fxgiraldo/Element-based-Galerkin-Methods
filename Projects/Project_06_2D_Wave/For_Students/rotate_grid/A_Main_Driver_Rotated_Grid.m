%---------------------------------------------------------------------%
%A sample driver for how to use the Grid Rotation routines for Project 4
%MA4245
%Written by F.X. Giraldo on 9/17/2014
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

%Input Data
nel=8; %Number of Elements
nop=4;    %Interpolation Order
noq=nop; %DO NOT CHANGE!
space_method='cg'; %CG or DG
plot_grid=0; %=0 don't plot, =1 plot
grid_rotation_angle=45; %CCW rotation in degrees

%Compute LGL Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);

%Compute Legendre Cardinal functions and derivatives
[psi,dpsi,xnq,wnq] = lagrange_basis(ngl,nq,xgl);

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

%Construct Initial Conditions and do Time-Integration
