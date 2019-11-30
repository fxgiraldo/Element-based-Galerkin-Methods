%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2010
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [mass, Jacobian] = create_mass_dg_fast(coord,nelem,ngl,nq,wnq,psi)

%Initialize
mass=zeros(ngl,ngl);
x=zeros(ngl,1);
Jacobian=zeros(nelem);

%Compute Jacobians
for e=1:nelem
   x(:)=coord(:,e);
   dx=x(ngl)-x(1);
   Jacobian(e)=dx/2;
end

%Compute 1 Mass Matrix
for l=1:nq
  wq=wnq(l);
  for i=1:ngl
     h_i=psi(i,l);
     for j=1:ngl
        h_j=psi(j,l);
        mass(i,j)=mass(i,j) + wq*h_i*h_j;
     end %j
  end %i
end %l



      