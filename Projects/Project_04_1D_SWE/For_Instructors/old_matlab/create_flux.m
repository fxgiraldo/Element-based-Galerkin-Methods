%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_flux(rhs,qp,qb,nelem,ngl,diss,gravity,delta_nl)

%Integrate Flux Terms
for e=1:nelem-1
   el=e;
   er=e+1;
   
   %LGL Integration
   hs_l=qp(1,ngl,el);
   U_l =qp(2,ngl,el);
   hb_l=qb(ngl,el);
   h_l =hs_l + hb_l;
   u_l =U_l/(h_l);
   
   hs_r=qp(1,1,er);
   U_r =qp(2,1,er);
   hb_r=qb(1,er);
   h_r =hs_r + hb_r;
   u_r =U_r/(h_r);
   
   %Compute Wave Speeds 
   clam_l=( abs(u_l) + sqrt(gravity*hs_l) )*delta_nl + sqrt(gravity*hb_l);
   clam_r=( abs(u_r) + sqrt(gravity*hs_r) )*delta_nl + sqrt(gravity*hb_r);
   clam=max(clam_l,clam_r);
   
   %Mass Flux
   f_l=U_l;
   f_r=U_r;
   flux_h=0.5*( f_l + f_r - diss*abs(clam)*(hs_r - hs_l) );
   
   %Momentum Flux
   f_l=( U_l^2/(h_l) + 0.5*gravity*hs_l^2 )*delta_nl + gravity*hs_l*hb_l;
   f_r=( U_r^2/(h_r) + 0.5*gravity*hs_r^2 )*delta_nl + gravity*hs_r*hb_r;
   flux_U=0.5*( f_l + f_r - diss*abs(clam)*(U_r - U_l));
   
   %Add to RHS to Left
   rhs(1,ngl,el)=rhs(1,ngl,el) - flux_h;
   rhs(2,ngl,el)=rhs(2,ngl,el) - flux_U;
   
   %Add to RHS to Right
   rhs(1,1,er)=rhs(1,1,er) + flux_h;
   rhs(2,1,er)=rhs(2,1,er) + flux_U;
end %ie

%Boundary Conditions

%Right Lateral Boundary
el=nelem;
hs_l=qp(1,ngl,el);
U_l =qp(2,ngl,el);
hb_l=qb(ngl,el);
h_l =hs_l + hb_l;
u_l =U_l/(h_l);

h_r =h_l;
hs_r=hs_l;
U_r =-U_l;
hb_r=hb_l;
u_r =U_r/(h_r);

%Compute Wave Speeds 
clam_l=( abs(u_l) + sqrt(gravity*hs_l) )*delta_nl + sqrt(gravity*hb_l);
clam_r=( abs(u_r) + sqrt(gravity*hs_r) )*delta_nl + sqrt(gravity*hb_r);
clam=max(clam_l,clam_r);

%Mass Flux
f_l=U_l;
f_r=U_r;
flux_h=0.5*( f_l + f_r - diss*abs(clam)*(hs_r - hs_l) );

%Momentum Flux
f_l=( U_l^2/(h_l) + 0.5*gravity*hs_l^2 )*delta_nl + gravity*hs_l*hb_l;
f_r=( U_r^2/(h_r) + 0.5*gravity*hs_r^2 )*delta_nl + gravity*hs_r*hb_r;
flux_U=0.5*( f_l + f_r - diss*abs(clam)*(U_r - U_l));

%Add to RHS of Left Element (Last element on Right Lateral Boundary)
rhs(1,ngl,el)=rhs(1,ngl,el) - flux_h;
rhs(2,ngl,el)=rhs(2,ngl,el) - flux_U;

%Left Lateral Boundary
er=1;
hs_r=qp(1,1,er);
U_r =qp(2,1,er);
hb_r=qb(1,er);
h_r =hs_r + hb_r;
u_r =U_r/(h_r);

h_l =h_r;
hs_l=hs_r;
U_l =-U_r;
hb_l=hb_r;
u_l =U_l/(h_l);

%Compute Wave Speeds 
clam_l=( abs(u_l) + sqrt(gravity*hs_l) )*delta_nl + sqrt(gravity*hb_l);
clam_r=( abs(u_r) + sqrt(gravity*hs_r) )*delta_nl + sqrt(gravity*hb_r);
clam=max(clam_l,clam_r);

%Mass Flux
f_l=U_l;
f_r=U_r;
flux_h=0.5*( f_l + f_r - diss*abs(clam)*(hs_r - hs_l) );

%Momentum Flux
f_l=( U_l^2/(h_l) + 0.5*gravity*hs_l^2 )*delta_nl + gravity*hs_l*hb_l;
f_r=( U_r^2/(h_r) + 0.5*gravity*hs_r^2 )*delta_nl + gravity*hs_r*hb_r;
flux_U=0.5*( f_l + f_r - diss*abs(clam)*(U_r - U_l));

%Add to RHS of Right Element (1st Element on Left Lateral Boundary)
rhs(1,1,er)=rhs(1,1,er) + flux_h;
rhs(2,1,er)=rhs(2,1,er) + flux_U;