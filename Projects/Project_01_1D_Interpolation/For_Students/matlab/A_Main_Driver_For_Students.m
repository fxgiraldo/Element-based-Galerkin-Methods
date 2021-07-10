%---------------------------------------------------------------------%
%This code computes the Interpolation using LGL points.
%Written by F.X. Giraldo on 4/2021
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

%Store Constants
N=4;
Npts=N+1;

%Compute Chebyshev Points
[xgl,wgl]=chebyshev_basis(Npts);
roots(:,1)=xgl;
weights(:,1)=wgl;
disp([' Chebyshev Points']);
xgl
wgl

%Compute Legendre Points
[xgl,wgl]=legendre_gauss(Npts);
roots(:,2)=xgl;
weights(:,2)=wgl;
disp([' Legendre Points']);
xgl
wgl

%Compute Lobatto Points
[xgl,wgl]=legendre_gauss_lobatto(Npts);
roots(:,3)=xgl;
weights(:,3)=wgl;
disp([' Lobatto Points']);
xgl
wgl

%Compute Equi-spaced Points
xgl=linspace(-1,1,Npts);
wgl=zeros(1,Npts);
roots(:,4)=xgl;
weights(:,4)=wgl;
disp([' Equi-spaced Points']);
xgl
wgl

plot(roots, weights,'o:','LineWidth',2)
legend('Chebyshev','Legendre','Lobatto','Equi-spaced')
set(gca, 'FontSize', 18);
