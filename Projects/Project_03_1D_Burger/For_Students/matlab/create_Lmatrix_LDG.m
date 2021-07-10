%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Lmatrix_LDG(rhs,Mmatrix,intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,q,alpha,beta,visc)

%Build Q: VOLUME
rhs_temp=create_Dmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,q);
%Build Q: FLUX
rhs_temp=create_Flux_q(rhs_temp,intma,nelem,ngl,q,alpha,beta);
%Build Q using Lift Operator
Q=Mmatrix\rhs_temp;

%Build q: VOLUME
rhs_temp=create_Dmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,Q);
%Build q: FLUX
rhs_temp=create_Flux_q(rhs_temp,intma,nelem,ngl,Q,beta,alpha);

%Add Viscous term
rhs = rhs + visc*rhs_temp;