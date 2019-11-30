%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function mass = create_mass_element(jac,wgl_matrix,nelem,ngl,element_order)

%Initialize
mass=zeros(ngl,nelem);

for e=1:nelem
    i=element_order(e);
    ii=i+1;
    mass(1:ii,e)=wgl_matrix(1:ii,i)*jac(e);
end %e
