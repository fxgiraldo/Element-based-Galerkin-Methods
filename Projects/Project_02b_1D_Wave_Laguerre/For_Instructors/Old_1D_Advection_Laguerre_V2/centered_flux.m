%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Naval Research Laboratory 
%           Monterey, CA 93943-5502
%---------------------------------------------------------------------%
function flux = centered_flux(q_l,q_r,diss)

flux=zeros(2,1);

%Left States
U_l=q_l(1);

%Right States
U_r=q_r(1);

%Compute Wave Speeds 
clam_l=abs(U_l);
clam_r=abs(U_r);
clam=max(clam_l,clam_r);
clam=0;

%Momentum Flux
f_l=0.5*U_l^2;
f_r=0.5*U_r^2;
flux(1)=0.5*( f_l + f_r - diss*abs(clam)*(U_r - U_l) );   
