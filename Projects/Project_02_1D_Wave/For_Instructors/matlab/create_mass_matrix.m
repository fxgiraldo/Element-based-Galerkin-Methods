%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on April 22, 2021
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function element_mass = create_mass_matrix(intma,coord,nelem,ngl,nq,wnq,psi)

%Initialize
element_mass=zeros(ngl,ngl,nelem);
x=zeros(ngl,1);

for e=1:nelem
   
   %Store Coordinates
   for i=1:ngl
       I=intma(i,e);
       x(i)=coord(I);
   end
   
   dx=x(ngl)-x(1);
   jac=dx/2;
      
   %Do LGL Integration
   for k=1:nq
      wk=wnq(k)*jac;
      for i=1:ngl
         h_i=psi(i,k);
         for j=1:ngl
            h_j=psi(j,k);
            element_mass(i,j,e)=element_mass(i,j,e) + wk*h_i*h_j;
         end %j
      end %i
   end %k
   
end %e



      