%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function mass = create_mass_dg_modal(nelem,ngl,nq,wnq,L,dx)

%Initialize
mass=zeros(ngl,ngl,nelem);

for ie=1:nelem
   jac=dx/2;
      
   %Do LGL Integration
   for l=1:nq
      wq=wnq(l)*jac;
      for i=1:ngl
         h_i=L(i,l);
         for j=1:ngl
            h_j=L(j,l);
            mass(i,j,ie)=mass(i,j,ie) + wq*h_i*h_j;
         end %j
      end %i
   end %l
end %ie



      
