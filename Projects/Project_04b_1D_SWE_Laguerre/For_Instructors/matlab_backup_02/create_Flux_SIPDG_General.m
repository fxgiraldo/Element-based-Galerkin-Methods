%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Flux_SIPDG_General(rhs,intma,coord,nelem,ngl,dpsi,q,tau,visc,lsponge)

%switches dPSIdx term in SIPDG
dd_L=0;
dd_R=0;
dd_int=0;

dd=dd_int;
% %FLUX
% for e=1:nelem-1
% 
%    %Left Values
%    dx_L=coord(ngl(e),e)-coord(1,e);
%    ksi_x_L=2/dx_L;
%    dhdx_L=dpsi(ngl(e),ngl(e),e)*ksi_x_L;
%    I=intma(ngl(e),e);
%    q_L=q(I,1);
%    U_L=q(I,2);
%    n_L=+1; %(pointing from Left -> Right)
%    h_visc_L=visc(I,1);
%    u_visc_L=visc(I,2);
%    dqdx_L=0;
%    dUdx_L=0;
%    for i=1:ngl(e)
%        I=intma(i,e);
%        dqdx_L=dqdx_L + dpsi(i,ngl(e),e)*q(I,1)*ksi_x_L;
%        dUdx_L=dUdx_L + dpsi(i,ngl(e),e)*q(I,2)*ksi_x_L;
%    end
%    q_L=q_L*h_visc_L;
%    U_L=U_L*u_visc_L;
%    dqdx_L=n_L*dqdx_L*h_visc_L;
%    dUdx_L=n_L*dUdx_L*u_visc_L;
%    dhqdx_L=n_L*dhdx_L*q_L;
%    dhUdx_L=n_L*dhdx_L*U_L;
% 
%    %Right Values
%    dx_R=coord(ngl(e),e+1)-coord(1,e+1);
%    ksi_x_R=2/dx_R;
%    dhdx_R=dpsi(1,1,e)*ksi_x_R;
%    I=intma(1,e+1);
%    q_R=q(I,1);
%    U_R=q(I,2);
%    n_R=-1;
%    h_visc_R=visc(I,1);
%    u_visc_R=visc(I,2);
%    dqdx_R=0;
%    dUdx_R=0;
%    for i=1:ngl(e)
%        I=intma(i,e+1);
%        dqdx_R=dqdx_R + dpsi(i,1,e)*q(I,1)*ksi_x_R;
%        dUdx_R=dUdx_R + dpsi(i,1,e)*q(I,2)*ksi_x_R;
%    end
%    q_R=q_R*h_visc_R;
%    U_R=U_R*u_visc_R;
%    dqdx_R=n_R*dqdx_R*h_visc_R;
%    dUdx_R=n_R*dUdx_R*u_visc_R;
%    dhqdx_R=n_R*dhdx_R*q_R;
%    dhUdx_R=n_R*dhdx_R*U_R;
% 
%    %Numerical Fluxes
%    mu_L=(ngl(e)-1)*(ngl(e)-1+1)/(dx_L)*tau; %Shabazi
%    mu_R=(ngl(e)-1)*(ngl(e)-1+1)/(dx_R)*tau;
%    mu=max(mu_L,mu_R);
%    dqdx_star=0.5*( dqdx_L - dqdx_R) - n_L*mu*(q_R-q_L);
%    dUdx_star=0.5*( dUdx_L - dUdx_R) - n_L*mu*(U_R-U_L);
%    dhqdx_star=0.5*(dhqdx_L - dhqdx_R);
%    dhUdx_star=0.5*(dhUdx_L - dhUdx_R);
% 
%    %Add to Left Element
%    I=intma(ngl(e),e);
%    rhs(I,1)=rhs(I,1) + dqdx_star + dd*(dhqdx_L-dhqdx_star);
%    rhs(I,2)=rhs(I,2) + dUdx_star + dd*(dhUdx_L-dhUdx_star);
% 
%    %Add to Right Element
%    I=intma(1,e+1);
%    rhs(I,1)=rhs(I,1) - dqdx_star + dd*(dhqdx_R-dhqdx_star);
%    rhs(I,2)=rhs(I,2) - dUdx_star + dd*(dhUdx_R-dhUdx_star);
% end %e=1,nelem-1

%------------------Left Boundary
dd=dd_L;

%Right Values (Interior)
e=1;
dx_R=coord(ngl(e),e)-coord(1,e);
ksi_x_R=2/dx_R;
dhdx_R=dpsi(1,1,e)*ksi_x_R;
I=intma(1,e);
q_R=q(I,1);
U_R=q(I,2);
n_R=-1; %(pointing from Right -> Left)
h_visc_R=visc(I,1);
u_visc_R=visc(I,2);
dqdx_R=0;
dUdx_R=0;
for i=1:ngl(e)
   I=intma(i,e);
   dqdx_R=dqdx_R + dpsi(i,1,e)*q(I,1)*ksi_x_R;
   dUdx_R=dUdx_R + dpsi(i,1,e)*q(I,2)*ksi_x_R;
end
q_R=q_R*h_visc_R;
U_R=U_R*u_visc_R;
dqdx_R=n_R*dqdx_R*h_visc_R;
dUdx_R=n_R*dUdx_R*u_visc_R;
dhqdx_R=n_R*dhdx_R*q_R;
dhUdx_R=n_R*dhdx_R*U_R;
   
%Left Values: No-Flux BC
n_L=+1;
dx_L=dx_R;
q_L=q_R;
U_L=-U_R;
dqdx_L=dqdx_R;
dUdx_L=dUdx_R;
dhqdx_L=-dhqdx_R;
dhUdx_L=dhUdx_R;

%Numerical Fluxes
mu_L=(ngl(e)-1)*(ngl(e)-1+1)/(dx_L)*tau; %Shabazi
mu_R=(ngl(e)-1)*(ngl(e)-1+1)/(dx_R)*tau;
mu=max(mu_L,mu_R);
dqdx_star=0.5*( dqdx_L - dqdx_R) - n_L*mu*(q_R-q_L);
dUdx_star=0.5*( dUdx_L - dUdx_R) - n_L*mu*(U_R-U_L);
dhqdx_star=0.5*(dhqdx_L - dhqdx_R);
dhUdx_star=0.5*(dhUdx_L - dhUdx_R);

%Add to Right Element
I=intma(1,e);
rhs(I,1)=rhs(I,1) - dqdx_star + dd*(dhqdx_R-dhqdx_star);
rhs(I,2)=rhs(I,2) - dUdx_star + dd*(dhUdx_R-dhUdx_star);
   
%------------------Right Boundary
if (lsponge==1)

    dd=dd_L;
    
    %Left Values (Interior)
    e=nelem;
    dx_L=coord(ngl(e),e)-coord(1,e);
    ksi_x_L=2/dx_L;
    dhdx_L=dpsi(ngl(e),ngl(e),e)*ksi_x_L;
    I=intma(ngl(e),e);
    q_L=q(I,1);
    U_L=q(I,2);
    n_L=+1; %(pointing from Left -> Right)
    h_visc_L=visc(I,1);
    u_visc_L=visc(I,2);
    dqdx_L=0;
    dUdx_L=0;
    for i=1:ngl(e)
        I=intma(i,e);
       dqdx_L=dqdx_L + dpsi(i,ngl(e),e)*q(I,1)*ksi_x_L;
       dUdx_L=dUdx_L + dpsi(i,ngl(e),e)*q(I,2)*ksi_x_L;
    end
    q_L=q_L*h_visc_L;
    U_L=U_L*u_visc_L;
    dqdx_L=n_L*dqdx_L*h_visc_L;
    dUdx_L=n_L*dUdx_L*u_visc_L;
    dhqdx_L=n_L*dhdx_L*q_L;
    dhUdx_L=n_L*dhdx_L*U_L;
    
    %Right Values: No-Flux BC
    q_R=q_L;
    U_R=-U_L;
    dqdx_R=dqdx_L;
    dUdx_R=dUdx_L;
    dhqdx_R=-dhqdx_L;
    dhUdx_R=dhUdx_L;
    
    %Numerical Fluxes
    mu_L=(ngl(e)-1)*(ngl(e)-1+1)/(dx_L)*tau; %Shabazi
    mu_R=(ngl(e)-1)*(ngl(e)-1+1)/(dx_R)*tau;
    mu=max(mu_L,mu_R);
    dqdx_star=0.5*( dqdx_L - dqdx_R) - n_L*mu*(q_R-q_L);
    dUdx_star=0.5*( dUdx_L - dUdx_R) - n_L*mu*(U_R-U_L);
    dhqdx_star=0.5*(dhqdx_L - dhqdx_R);
    dhUdx_star=0.5*(dhUdx_L - dhUdx_R);
    
    %Add to Left Element
    I=intma(ngl(e),e);
    rhs(I,1)=rhs(I,1) + dqdx_star + dd*(dhqdx_L-dhqdx_star);
    rhs(I,2)=rhs(I,2) + dUdx_star + dd*(dhUdx_L-dhUdx_star);

end
