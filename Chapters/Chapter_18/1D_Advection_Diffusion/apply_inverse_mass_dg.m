%---------------------------------------------------------------------%
%This function computes the Mass Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_inverse_mass_dg(rhs,mass_inv,nelem,ngl)

%Multiply by Inverse Mass matrix
for e=1:nelem
 for j=1:ngl
    for i=1:ngl
       mass_t(i,j)=mass_inv(i,j,e);
    end %i
 end %j
 rhs_t=rhs(:,e);
 rhs(:,e)=mass_t*rhs_t;
end %ie
      
