%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Lmatrix_LDG(Mmatrix,intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,q,qe,qe_x)

%Initialize
rhs=zeros(npoin,1);
Q=zeros(npoin,1);
alpha=0.5; beta=1.0-alpha;

%Build Q: VOLUME
rhs=create_Dmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,q);
%Build Q: FLUX
rhs=create_Flux_q(rhs,intma,nelem,ngl,q,qe,alpha,beta);
%Build Q using Lift Operator
Q=Mmatrix\rhs;

%Build q: VOLUME
rhs=create_Dmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,Q);
%Build q: FLUX
rhs=create_Flux_q(rhs,intma,nelem,ngl,Q,qe_x,beta,alpha);