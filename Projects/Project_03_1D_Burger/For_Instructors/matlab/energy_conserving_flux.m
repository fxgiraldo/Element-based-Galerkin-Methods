%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Naval Research Laboratory 
%           Monterey, CA 93943-5502
%---------------------------------------------------------------------%
function flux = energy_conserving_flux(q_l,q_r,diss)

flux=zeros(2,1);
eps=1e-1;

%Left States
U_l=q_l(1);

%Right States
U_r=q_r(1);

%Compute Wave Speeds 
clam_l=abs(U_l);
clam_r=abs(U_r);
clam=0.5*max(clam_l,clam_r); %Rusanov
clam=0.5*abs(U_l + U_r); %Roe
clam=1.0/12.0*(U_r - U_l); %Energy

%Momentum Flux
f_l=0.5*U_l^2;
f_r=0.5*U_r^2;
flux(1)=0.5*( f_l + f_r) - diss*clam*(U_r - U_l);   


%Momentum Flux
% f_l=0.5*U_l^2;
% f_r=0.5*U_r^2;
% flux(1)=0.5*( f_l + f_r - diss*1.0/6.0*(U_r - U_l)^2 ); 
% %flux(1)=1.0/6.0*( U_l^2 + U_l*U_r + U_r^2 ); 
