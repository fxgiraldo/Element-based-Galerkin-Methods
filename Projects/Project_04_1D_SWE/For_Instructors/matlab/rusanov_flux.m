%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function flux = rusanov_flux(q_l,q_r,diss,gravity,delta_nl)

flux=zeros(2,1);

%Left States
hs_l=q_l(1);
hb_l=q_l(2);
u_l=q_l(3);
h_l=hs_l + hb_l;
U_l=h_l*u_l;

%Right States
hs_r=q_r(1);
hb_r=q_r(2);
u_r=q_r(3);
h_r=hs_r + hb_r;
U_r=h_r*u_r;

%Compute Wave Speeds 
clam_l=( abs(u_l) + sqrt(gravity*hs_l) )*delta_nl + sqrt(gravity*hb_l);
clam_r=( abs(u_r) + sqrt(gravity*hs_r) )*delta_nl + sqrt(gravity*hb_r);
clam=max(clam_l,clam_r);

%Mass Flux
f_l=U_l;
f_r=U_r;
flux(1)=0.5*( f_l + f_r - diss*abs(clam)*(hs_r - hs_l) );

%Momentum Flux
f_l=( U_l^2/(h_l) + 0.5*gravity*hs_l^2 )*delta_nl + gravity*hs_l*hb_l;
f_r=( U_r^2/(h_r) + 0.5*gravity*hs_r^2 )*delta_nl + gravity*hs_r*hb_r;
flux(2)=0.5*( f_l + f_r - diss*abs(clam)*(U_r - U_l));
   
