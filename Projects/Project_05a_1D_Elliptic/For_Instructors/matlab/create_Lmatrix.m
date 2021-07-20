%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Lmatrix = create_Lmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,iperiodic)

%Initialize
Lmatrix=zeros(npoin,npoin);

for e=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      ip=intma(e,i);
      inode(i)=iperiodic(ip);
      x(i)=coord(ip);
   end
   
   dx=x(ngl)-x(1);
   jac=dx/2;
   dksi_dx=2/dx;
   
   %LGL Integration
   for i=1:ngl
       ip=inode(i);
       for j=1:ngl
           jp=inode(j);
           for k=1:nq
               wq=wnq(k)*jac;
               dhdx_ik=dpsi(i,k)*dksi_dx;
               dhdx_jk=dpsi(j,k)*dksi_dx;
               Lmatrix(ip,jp)=Lmatrix(ip,jp) + wq*dhdx_ik*dhdx_jk;
           end %k
       end %j
   end %i
end %e

%Impose Dirichlet Boundary Conditions
Lmatrix(1,:)=0;
Lmatrix(1,1)=1;
Lmatrix(npoin,:)=0;
Lmatrix(npoin,npoin)=1;