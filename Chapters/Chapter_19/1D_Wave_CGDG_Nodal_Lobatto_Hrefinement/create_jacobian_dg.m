%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function jac = create_jacobian_dg(coord,nelem,ngl)

%Initialize
jac=zeros(nelem,1);

%Store Element Jacobian
for e=1:nelem
   dx=coord(ngl,e)-coord(1,e);
   jac(e)=dx/2;
end %e
