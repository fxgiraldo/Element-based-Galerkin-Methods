%---------------------------------------------------------------------%
%This function computes the 2D Mass Matrix on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Mmatrix = create_Mmatrix(intma,jac,wnq,psi,iperiodic,npoin,...
                   nelem,ngl,nq)

%Initialize
Mmatrix=zeros(npoin,npoin);
inode=zeros(ngl,ngl);

for e=1:nelem
   
   %Store Coordinates
   for j=1:ngl
   for i=1:ngl
      ip=intma(i,j,e);
      inode(i,j)=iperiodic(ip);
   end %i
   end %j
   
   %Do LGL Integration
   for l=1:nq
   for k=1:nq
       wq=wnq(k)*wnq(l)*jac(k,l,e);
       
       %Loop through I points
       for j=1:ngl
       for i=1:ngl
           ip=inode(i,j);
           h_i=psi(i,k)*psi(j,l);
           
           %Loop through J points
           for n=1:ngl
           for m=1:ngl
               jp=inode(m,n);
               h_j=psi(m,k)*psi(n,l);
               Mmatrix(ip,jp)=Mmatrix(ip,jp) + wq*h_i*h_j;
           end %m
           end %n
       end %i
       end%j
   end %k
   end %l
end %e

%Periodicity
for i=1:npoin
    j=iperiodic(i);
    if (i ~= j) 
        Mmatrix(i,:)=0;
        Mmatrix(i,i)=1;
    end
end


      
