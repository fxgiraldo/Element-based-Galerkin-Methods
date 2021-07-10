%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Flux_SIPDG(rhs,intma,coord,nelem,ngl,nq,psi,dpsi,q,tau)

%FLUX
for e=1:nelem
    
   eL=e;
   eR=e+1;
   if (eR > nelem)
       eR=1;
   end
   
   %Left Values
   dx_L=coord(ngl,eL)-coord(1,eL);
   jac_L=dx_L/2;
   ksi_x_L=2/dx_L;
   IL=intma(ngl,eL);
   dhdx_L=dpsi(ngl,nq)*ksi_x_L;
   h_L=psi(ngl,nq);
   q_L=q(IL);
   dqdx_L=0;
   for i=1:ngl
       ip=intma(i,eL);
       dqdx_L=dqdx_L + dpsi(i,nq)*q(ip)*ksi_x_L;
   end
   n_L=+1; %(pointing from Left -> Right)
   
   %Right Values
   dx_R=coord(ngl,eR)-coord(1,eR);
   ksi_x_R=2/dx_R;
   IR=intma(1,eR);
   dhdx_R=dpsi(1,1)*ksi_x_R;
   h_R=psi(1,1);
   q_R=q(IR);
   dqdx_R=0;
   for i=1:ngl
       ip=intma(i,eR);
       dqdx_R=dqdx_R + dpsi(i,1)*q(ip)*ksi_x_R;
   end
   n_R=-1;
   
   %Numerical Fluxes
   q_star=0.5*(q_L + q_R);
   mu_L=(ngl-1)*(ngl-1+1)/(dx_L)*tau; %Shabazi
   mu_R=(ngl-1)*(ngl-1+1)/(dx_R)*tau;
   mu=max(mu_L,mu_R);
   dqdx_mean=0.5*( dqdx_L + dqdx_R); 
   dqdx_diss=mu*(q_R-q_L);
   
   %Add to Left/Right Elements
   rhs(IL)=rhs(IL) + n_L*h_L*dqdx_mean - 0*h_L*dqdx_diss + 0*n_L*dhdx_L*(q_L-q_star); %not valid for periodic
   rhs(IR)=rhs(IR) + n_R*h_R*dqdx_mean - 0*h_R*dqdx_diss + 0*n_R*dhdx_R*(q_R-q_star);
end