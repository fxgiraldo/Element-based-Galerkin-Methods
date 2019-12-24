%----------------------------------------------------------------------%
%This subroutine interpolates edge values onto Quadrature Points and builds 
%derivatives of the function at the edges 
%Written by Francis X. Giraldo on 8/2015
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function [q_k, q_x, q_y] = compute_flux_SIP_qterms_Inexact(q,psi,dpsi,ngl,nq,...
                           ksi_x,ksi_y,eta_x,eta_y,iloc,e,l,imap)

%Interpolate Left-State onto Quadrature Points
q_k=0;q_e=0;q_n=0;q_x=0; q_y=0;

q_k=q(l);
kl=imap(iloc,1,l);
ll=imap(iloc,2,l);

for i=1:ngl
  il=imap(iloc,1,i);
  jl=imap(iloc,2,i);
  h_j=psi(il,kl);
  h_k=psi(jl,ll);
  h_jk=h_j*h_k;
  dhde_jk=dpsi(il,kl)*psi(jl,ll);
  dhdn_jk=psi(il,kl)*dpsi(jl,ll);
  q_e=q_e + dhde_jk*q(i);
  q_n=q_n + dhdn_jk*q(i);
end %i 
q_x=q_e*ksi_x(e,kl,ll) + q_n*eta_x(e,kl,ll);
q_y=q_e*ksi_y(e,kl,ll) + q_n*eta_y(e,kl,ll);
