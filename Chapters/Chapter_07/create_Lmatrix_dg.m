%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function laplacian_matrix = create_Lmatrix_dg(coord,nelem,ngl,nq,wnq,dpsi)

%Initialize
laplacian_matrix=zeros(ngl,ngl,nelem);

for e=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      x(i)=coord(i,e);
   end
   
   dx=x(ngl)-x(1);
   jac=dx/2;
   dksi_dx=2/dx;
   
   %LGL Integration
   for i=1:ngl
       for j=1:ngl
           for k=1:nq
               wq=wnq(k)*jac;
               dhdx_ik=dpsi(i,k)*dksi_dx;
               dhdx_jk=dpsi(j,k)*dksi_dx;
               laplacian_matrix(i,j,e)=laplacian_matrix(i,j,e) + wq*dhdx_ik*dhdx_jk;
           end %k
       end %j
   end %i
end %e
