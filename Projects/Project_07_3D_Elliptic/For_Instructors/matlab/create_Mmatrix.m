%---------------------------------------------------------------------%
%This function computes the 2D Mass Matrix on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Mmatrix = create_Mmatrix(intma,jac,wnq,psi,npoin,nelem,ngl,nq)

%Initialize
Mmatrix=zeros(npoin,npoin);
inode=zeros(ngl,ngl,ngl);

for e=1:nelem
   
   %Store Coordinates
   for k=1:ngl
   for j=1:ngl
   for i=1:ngl
      I=intma(i,j,k,e);
      inode(i,j,k)=I;
   end %i
   end %j
   end %k
   
   %Do LGL Integration
   for n=1:nq
   for m=1:nq
   for l=1:nq
       wq=wnq(l)*wnq(m)*wnq(n)*jac(l,m,n,e);
       
       %Loop through I points
       for k=1:ngl
       for j=1:ngl
       for i=1:ngl
           I=inode(i,j,k);
           h_i=psi(i,l)*psi(j,m)*psi(k,n);
           
           %Loop through J points
           for kj=1:ngl
           for jj=1:ngl
           for ij=1:ngl
               J=inode(ij,jj,kj);
               h_j=psi(ij,l)*psi(jj,m)*psi(kj,n);
               Mmatrix(I,J)=Mmatrix(I,J) + wq*h_i*h_j;
           end %ij
           end %jj
           end %kj
       end %i
       end %j
       end %k
   end %l
   end %m
   end %n
end %e