%----------------------------------------------------------------------%
%This subroutine interpolates edge values onto Quadrature Points and builds 
%derivatives of the function at the edges 
%Written by Francis X. Giraldo on 8/2015
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function [q_k] = compute_flux_LDG_qterms(q,psi,ngl,l)

%Interpolate Left-State onto Quadrature Points
q_k=0;
for i=1:ngl
    q_k=q_k + psi(i,l)*q(i);
end %i 

        
