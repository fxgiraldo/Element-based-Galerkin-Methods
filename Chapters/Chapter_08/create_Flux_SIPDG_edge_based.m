%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Flux_SIPDG_edge_based(rhs,intma,coord,nelem,ngl,nq,psi,dpsi,q,qe,qe_x,tau,sipdg_flag)

%FLUX
for e=1:nelem-1
    
   %Left Values
   dx_L=coord(ngl,e)-coord(1,e);
   jac_L=dx_L/2;
   ksi_x_L=2/dx_L;
   IL=intma(ngl,e);
   dhdx_L=dpsi(ngl,nq)*ksi_x_L;
   h_L=psi(ngl,nq);
   q_L=q(IL);
   dqdx_L=0;
   for i=1:ngl
       ip=intma(i,e);
       dqdx_L=dqdx_L + dpsi(i,nq)*q(ip)*ksi_x_L;
   end
   n_L=+1; %(pointing from Left -> Right)
   
   %Right Values
   dx_R=coord(ngl,e+1)-coord(1,e+1);
   ksi_x_R=2/dx_R;
   IR=intma(1,e+1);
   dhdx_R=dpsi(1,1)*ksi_x_R;
   h_R=psi(1,1);
   q_R=q(IR);
   dqdx_R=0;
   for i=1:ngl
       ip=intma(i,e+1);
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
   
   %Add to Left Element
   rhs(IL)=rhs(IL) + n_L*h_L*dqdx_mean - h_L*dqdx_diss + n_L*dhdx_L*(q_L-q_star);
   rhs(IR)=rhs(IR) + n_R*h_R*dqdx_mean - h_R*dqdx_diss + n_R*dhdx_R*(q_R-q_star);
end

%------------------Left Boundary
%Left Values
dx=coord(ngl,1)-coord(1,1);
jac=dx/2;
ksi_x=2/dx;
IL=intma(1,1);
dhdx_L=dpsi(1,1)*ksi_x;
h_L=psi(1,1);
q_L=q(IL);
dqdx_L=0;
for i=1:ngl
   ip=intma(i,1);
   dqdx_L=dqdx_L + dpsi(i,1)*q(ip)*ksi_x;
end
n_L=-1; %(pointing from Left -> Right)

%Right Values (Boundary)
q_R=qe(intma(1,1));
dqdx_R=qe_x(intma(1,1));

%Numerical Fluxes
q_star=0.5*(q_L + q_R);
mu=(ngl-1)*(ngl-1+1)/(dx)*tau; %Shabazi
dqdx_star=0.5*( dqdx_L + dqdx_R) - mu*(q_R-q_L);

%Add to Left Element
rhs(IL)=rhs(IL) + n_L*h_L*dqdx_star + n_L*dhdx_L*(q_L-q_star);

%------------------Right Boundary
%Left Values
dx=coord(ngl,nelem)-coord(1,nelem);
jac=dx/2;
ksi_x=2/dx;
IL=intma(ngl,nelem);
dhdx_L=dpsi(ngl,nq)*ksi_x;
h_L=psi(ngl,nq);
q_L=q(IL);
dqdx_L=0;
for i=1:ngl
   ip=intma(i,nelem);
   dqdx_L=dqdx_L + dpsi(i,nq)*q(ip)*ksi_x;
end
n_L=+1; %(pointing from Left -> Right)

%Right Values
q_R=qe(intma(ngl,nelem));
dqdx_R=qe_x(intma(ngl,nelem));

%Numerical Fluxes
q_star=0.5*(q_L + q_R);
mu=(ngl-1)*(ngl-1+1)/(dx)*tau; %Shabazi
dqdx_star=0.5*( dqdx_L + dqdx_R) - mu*(q_R-q_L);

%Add to Element
rhs(IL)=rhs(IL) + n_L*h_L*dqdx_star + n_L*dhdx_L*(q_L-q_star);
