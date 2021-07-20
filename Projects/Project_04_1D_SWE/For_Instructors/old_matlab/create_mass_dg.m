%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Naval Research Laboratory 
%           Monterey, CA 93943-5502
%---------------------------------------------------------------------%
function mass = create_mass_dg(coord,nelem,ngl,nq,wnq,psi)

%Initialize
mass=zeros(ngl,ngl,nelem);

for ie=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      x(i)=coord(i,ie);
   end
   
   dx=x(ngl)-x(1);
   jac=dx/2;
      
   %Do LGL Integration
   for l=1:nq
      wq=wnq(l)*jac;
      for i=1:ngl
         h_i=psi(i,l);
         for j=1:ngl
            h_j=psi(j,l);
            mass(i,j,ie)=mass(i,j,ie) + wq*h_i*h_j;
         end %j
      end %i
   end %l
end %ie



      