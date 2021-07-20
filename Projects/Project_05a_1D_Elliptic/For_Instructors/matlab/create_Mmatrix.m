%---------------------------------------------------------------------%
%This function computes the Mass Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Mmatrix = create_Mmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,iperiodic)

%Initialize
Mmatrix=zeros(npoin,npoin);

for ie=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      ip=intma(ie,i);
      inode(i)=iperiodic(ip);
      x(i)=coord(ip);
   end
   
   dx=x(ngl)-x(1);
   jac=dx/2;
   
   %Do LGL Integration
   for i=1:ngl
      ip=inode(i);
      for j=1:ngl
         jp=inode(j);
         for k=1:nq
            wq=wnq(k)*jac;
            h_ik=psi(i,k);
            h_jk=psi(j,k);
            Mmatrix(ip,jp)=Mmatrix(ip,jp) + wq*h_ik*h_jk;
         end %k
      end %j
   end %i
end %ie

%Impose Dirichlet Boundary Conditions
Mmatrix(1,:)=0; Mmatrix(1,1)=1;
Mmatrix(npoin,:)=0; Mmatrix(npoin,npoin)=1;

      
