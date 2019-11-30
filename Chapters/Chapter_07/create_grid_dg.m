%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [coord,intma] = create_grid_dg(ngl,nelem,xgl)

%Set some constants
% xmin=-1;
xmin=-1;
xmax=+1;
dx=(xmax-xmin)/nelem;

%Generate Grid Points
ip=1;
coord(1,1)=xmin;
for i=1:nelem
   x0=xmin + (i-1)*dx;
   coord(1,i)=x0;
   intma(1,i)=ip;
   for j=2:ngl
      ip=ip + 1;
      coord(j,i)=( xgl(j)+1 )*dx/2 + x0;
      intma(j,i)=ip;
   end %j
end %i
      


      