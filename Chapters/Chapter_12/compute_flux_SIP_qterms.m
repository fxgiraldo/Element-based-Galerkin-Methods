%----------------------------------------------------------------------%
%This subroutine interpolates edge values onto Quadrature Points and builds 
%derivatives of the function at the edges 
%Written by Francis X. Giraldo on 8/2015
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function [q_k, q_x, q_y] = compute_flux_SIP_qterms(q,psi,dpsi,ngl,nq,...
                           ksi_x,ksi_y,eta_x,eta_y,iloc,e,l)

%Interpolate Left-State onto Quadrature Points
q_k=0;q_e=0;q_n=0;
switch (iloc)
case (1)
%Side 1: eta=-1
    for i=1:ngl
      h_j=psi(i,l);
      h_k=psi(1,1);
      h_jk=h_j*h_k;
      dhde_jk=dpsi(i,l)*psi(1,1);
      dhdn_jk=psi(i,l)*dpsi(1,1);
      q_k=q_k + h_jk*q(i);
      q_e=q_e + dhde_jk*q(i);
      q_n=q_n + dhdn_jk*q(i);
    end %i 
    q_x=q_e*ksi_x(e,l,1) + q_n*eta_x(e,l,1);
    q_y=q_e*ksi_y(e,l,1) + q_n*eta_y(e,l,1);
case (2)
%Side 2: ksi=+1
    for i=1:ngl
      h_j=psi(ngl,nq);
      h_k=psi(i,l);
      h_jk=h_j*h_k;
      dhde_jk=dpsi(ngl,nq)*psi(i,l);
      dhdn_jk=psi(ngl,nq)*dpsi(i,l);
      q_k=q_k + h_jk*q(i);
      q_e=q_e + dhde_jk*q(i);
      q_n=q_n + dhdn_jk*q(i);
    end %i 
    q_x=q_e*ksi_x(e,nq,l) + q_n*eta_x(e,nq,l);
    q_y=q_e*ksi_y(e,nq,l) + q_n*eta_y(e,nq,l);
case (3)
    %Side 3: eta=+1
    for i=1:ngl
      h_j=psi(i,l);
      h_k=psi(ngl,nq);
      h_jk=h_j*h_k;
      dhde_jk=dpsi(i,l)*psi(ngl,nq);
      dhdn_jk=psi(i,l)*dpsi(ngl,nq);
      q_k=q_k + h_jk*q(i);
      q_e=q_e + dhde_jk*q(i);
      q_n=q_n + dhdn_jk*q(i);
    end %i  
    q_x=q_e*ksi_x(e,l,nq) + q_n*eta_x(e,l,nq);
    q_y=q_e*ksi_y(e,l,nq) + q_n*eta_y(e,l,nq);
case (4)
    %Side 4: ksi=-1
    for i=1:ngl
      h_j=psi(1,1);
      h_k=psi(i,l);
      h_jk=h_j*h_k;
      dhde_jk=dpsi(1,1)*psi(i,l);
      dhdn_jk=psi(1,1)*dpsi(i,l);
      q_k=q_k + h_jk*q(i);
      q_e=q_e + dhde_jk*q(i);
      q_n=q_n + dhdn_jk*q(i);
    end %i 
    q_x=q_e*ksi_x(e,1,l) + q_n*eta_x(e,1,l);
    q_y=q_e*ksi_y(e,1,l) + q_n*eta_y(e,1,l);
end %switch 

        
