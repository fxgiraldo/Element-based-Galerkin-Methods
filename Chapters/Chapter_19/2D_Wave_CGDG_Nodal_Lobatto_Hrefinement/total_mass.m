%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function mass = total_mass(q,jac,sfc,nsfc,ngl)

mass = 0;
for is=1:nsfc
   ie = sfc(is);
   for j=1:ngl
       for i=1:ngl
         mass = mass + jac(ie,i,j)*q(ie,i,j);
       end
   end
end

end
