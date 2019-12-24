%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix using
%the Symmetric Interior Penalty Method
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Lmatrix_SIPDG(rhs,intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,q,tau,visc)

rhs_temp=zeros(npoin,1);
%VOLUME Term: IBP
rhs_temp=create_Lmatrix_IBP(rhs_temp,intma,coord,nelem,ngl,nq,wnq,dpsi,q);
%Flux
rhs_temp=create_Flux_SIPDG(rhs_temp,intma,coord,nelem,ngl,nq,psi,dpsi,q,tau);
%Update
rhs=rhs + visc*rhs_temp;