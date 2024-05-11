%---------------------------------------------------------------------%
%This function computes the Derivative Mapping for General 2D Quad Grids.
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [f_ksi,f_eta] = map_deriv(psi,dpsi,f,ngl,nq)

for l=1:nq
for k=1:nq
   
    sum_ksi=0;
    sum_eta=0;

    for j=1:ngl
    for i=1:ngl
        sum_ksi=sum_ksi + dpsi(i,k)*psi(j,l)*f(i,j);
        sum_eta=sum_eta + psi(i,k)*dpsi(j,l)*f(i,j);
    end %i
    end %j
   
    f_ksi(k,l)=sum_ksi;
    f_eta(k,l)=sum_eta;

end %k
end %l
