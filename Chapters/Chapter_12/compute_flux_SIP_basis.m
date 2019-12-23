%----------------------------------------------------------------------%
%This subroutine interpolates basis functions onto Quadrature Points 
%Written by Francis X. Giraldo on 8/2015
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function [psi_x, psi_y] = compute_flux_SIP_basis(psi,dpsi,ngl,nq,...
                           ksi_x,ksi_y,eta_x,eta_y,iloc,e,i,l)

%Interpolate State onto Quadrature Points
switch (iloc)
case (1)
%Side 1: eta=-1
    psi_e=dpsi(i,l)*psi(1,1);
    psi_n=psi(i,l)*dpsi(1,1);
    psi_x=psi_e*ksi_x(e,l,1) + psi_n*eta_x(e,l,1);
    psi_y=psi_e*ksi_y(e,l,1) + psi_n*eta_y(e,l,1);
case (2)
%Side 2: ksi=+1
    psi_e=dpsi(ngl,nq)*psi(i,l);
    psi_n=psi(ngl,nq)*dpsi(i,l);
    psi_x=psi_e*ksi_x(e,nq,l) + psi_n*eta_x(e,nq,l);
    psi_y=psi_e*ksi_y(e,nq,l) + psi_n*eta_y(e,nq,l);
case (3)
    %Side 3: eta=+1
    psi_e=dpsi(i,l)*psi(ngl,nq);
    psi_n=psi(i,l)*dpsi(ngl,nq);
    psi_x=psi_e*ksi_x(e,l,nq) + psi_n*eta_x(e,l,nq);
    psi_y=psi_e*ksi_y(e,l,nq) + psi_n*eta_y(e,l,nq);
case (4)
    %Side 4: ksi=-1
    psi_e=dpsi(1,1)*psi(i,l);
    psi_n=psi(1,1)*dpsi(i,l);
    psi_x=psi_e*ksi_x(e,1,l) + psi_n*eta_x(e,1,l);
    psi_y=psi_e*ksi_y(e,1,l) + psi_n*eta_y(e,1,l);
end %switch 

        