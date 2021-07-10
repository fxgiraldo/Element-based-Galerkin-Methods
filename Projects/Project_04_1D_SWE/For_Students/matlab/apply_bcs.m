%---------------------------------------------------------------------%
%This function computes the Mass Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_bcs(rhs,qp,qb,gravity,nelem,ngl,delta_nl)
%No-Flux BC
%Boundary Integral

%Right Lateral Boundary
ps_k=qp(1,ngl,nelem);
b_k =qb(ngl,nelem);
p_k=ps_k + b_k;
pu_k=qp(2,ngl,nelem);
flux_p=0*pu_k;
flux_pu=( 0*pu_k^2/(p_k) + 0.5*gravity*ps_k^2 )*delta_nl + gravity*b_k*ps_k;
rhs(1,ngl,nelem)=rhs(1,ngl,nelem) - flux_p;
rhs(2,ngl,nelem)=rhs(2,ngl,nelem) - flux_pu;

%Left Lateral Boundary
ps_k=qp(1,1,1);
b_k =qb(1,1);
pu_k=qp(2,1,1);
flux_p=0*pu_k;
flux_pu=( 0*pu_k^2/(p_k) + 0.5*gravity*ps_k^2 )*delta_nl + gravity*b_k*ps_k;
rhs(1,1,1)=rhs(1,1,1) + flux_p;
rhs(2,1,1)=rhs(2,1,1) + flux_pu;





      
