%---------------------------------------------------------------------%
%This function computes the Mass Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_communicator(rhs,space_method,intma,intma_cg,coord,Mmatrix,npoin,npoin_cg,nelem,ngl,nq,wnq,psi,dpsi,alpha,ivolume_integrate)

%Apply Communicator
if ( strcmp(space_method,'cgd') > 0 ) %CGd
  rhs = apply_dss(rhs,intma,intma_cg,coord,Mmatrix,npoin,npoin_cg,nelem,ngl,nq,wnq,psi,dpsi,alpha,ivolume_integrate);
else %CGC or DG         
  %Multiply by Inverse Mass matrix
  rhs(:,1)=Mmatrix\rhs(:,1);
  rhs(:,2)=Mmatrix\rhs(:,2);
end