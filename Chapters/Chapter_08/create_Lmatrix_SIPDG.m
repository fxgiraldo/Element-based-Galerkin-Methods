%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix using
%the Symmetric Interior Penalty Method
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Lmatrix_SIPDG(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,q,qe,qe_x,tau,sipdg_flag)

%Initialize
rhs=zeros(npoin,1);

%VOLUME Term: IBP
rhs=create_Lmatrix_IBP(intma,coord,npoin,nelem,ngl,nq,wnq,dpsi,q);
%Flux
%rhs=create_Flux_SIPDG(rhs,intma,coord,nelem,ngl,nq,psi,dpsi,q,qe,qe_x,tau,sipdg_flag);
rhs=create_Flux_SIPDG_edge_based(rhs,intma,coord,nelem,ngl,nq,psi,dpsi,q,qe,qe_x,tau,sipdg_flag);
