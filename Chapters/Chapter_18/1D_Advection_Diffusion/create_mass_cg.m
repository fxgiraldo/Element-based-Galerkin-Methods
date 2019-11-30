%---------------------------------------------------------------------%
%This function computes the Mass Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Mmatrix = create_mass_cg(intma,coord,npoin,nelem,ngl,nq,wnq,psi,iperiodic)

%Initialize
Mmatrix=zeros(npoin,npoin);
inode=zeros(ngl,1);

for ie=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      ip=intma(i,ie);
      inode(i)=iperiodic(ip);
      x(i)=coord(i,ie);
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
            h_i=psi(i,k);
            h_j=psi(j,k);
            Mmatrix(ip,jp)=Mmatrix(ip,jp) + wq*h_i*h_j;
         end %k
      end %j
   end %i
end %ie

%Periodicity
Mmatrix(npoin,npoin)=1;

      
