%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Dhatmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,q,qe,alpha,beta)

%Initialize
rhs=zeros(npoin,1);

%Build Q: VOLUME
rhs=create_Dmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,q);
%Build Q: FLUX
% rhs=create_Flux_LDG(rhs,intma,nelem,ngl,q,qe,alpha,beta);
rhs=create_Flux_LDG_element(rhs,intma,nelem,ngl,q,qe,alpha,beta);