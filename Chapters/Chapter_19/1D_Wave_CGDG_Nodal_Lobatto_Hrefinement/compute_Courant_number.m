%---------------------------------------------------------------------%
%This function computes the Courant number.
%Written by F.X. Giraldo on 10/2003
%           Naval Postgraduate School
%           Monterey, CA 93943
%---------------------------------------------------------------------%
function [Courant,dx] = compute_Courant_number(coord,ngl,u,sfc,nsfc,dt)

%Find Minimum spacing
dx=coord(ngl,1)-coord(1,1);

for ee=1:nsfc
    e=sfc(ee);
    x0=coord(1,e);
    xN=coord(2,e);
    dx=min(xN-x0,dx);
end %e

Courant=u*dt/dx;
