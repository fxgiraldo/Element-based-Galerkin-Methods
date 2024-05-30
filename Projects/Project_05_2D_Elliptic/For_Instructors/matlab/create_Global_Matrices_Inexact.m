%---------------------------------------------------------------------%
%This function computes the 2D Laplacian Matrix on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%Optimized by F.X. Giraldo on May 26, 2024 to reduce cost to O(N^2d)
%instead of O(N^3d)
%---------------------------------------------------------------------%
function [Mmatrix,Lmatrix] = create_Global_Matrices_Inexact(intma,jac,wnq,ksi_x,ksi_y,eta_x,...         
                   eta_y,psi,dpsi,iperiodic,npoin,nelem,ngl,nq)

%Initialize
Mmatrix=zeros(npoin,npoin);
Lmatrix=zeros(npoin,npoin);
inode=zeros(ngl,ngl);

for e=1:nelem
   
   %Store Coordinates
   for j=1:ngl
   for i=1:ngl
      I=intma(i,j,e);
      inode(i,j)=iperiodic(I);
   end %i
   end %j
   
   %Do LGL Integration
   for l=1:nq
   for k=1:nq
       wq=wnq(k)*wnq(l)*jac(k,l,e);
       e_x=ksi_x(k,l,e);
       e_y=ksi_y(k,l,e);
       n_x=eta_x(k,l,e);
       n_y=eta_y(k,l,e);

       %Mass Matrix
       I=inode(k,l);
       Mmatrix(I,I)=Mmatrix(I,I) + wq;

       %Laplacian Matrix
       %Loop through I points
       for i=1:ngl

           %dXI derivatives
           I=inode(i,l); %j=l
           h_e=dpsi(i,k)*psi(l,l); %j=l
           dhdx_i=h_e*e_x ;
           dhdy_i=h_e*e_y ;
             
           for m=1:ngl
               J=inode(m,l); %n=l
               h_e=dpsi(m,k)*psi(l,l); %n=l
               dhdx_j=h_e*e_x;
               dhdy_j=h_e*e_y;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j);
               J=inode(k,m); %m=k and swap n-> m
               h_n=psi(k,k)*dpsi(m,l);
               dhdx_j=h_n*n_x;
               dhdy_j=h_n*n_y;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j);
           end %m

           %dETA derivatives
           I=inode(k,i); %i=k and swapped j->i
           h_n=psi(k,k)*dpsi(i,l); %i=k and swapped j->i
           dhdx_i=h_n*n_x;
           dhdy_i=h_n*n_y;
             
           for m=1:ngl
               J=inode(m,l); %n=l
               h_e=dpsi(m,k)*psi(l,l); %n=l
               dhdx_j=h_e*e_x;
               dhdy_j=h_e*e_y;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j);
               J=inode(k,m); %m=k and swap n-> m
               h_n=psi(k,k)*dpsi(m,l);
               dhdx_j=h_n*n_x;
               dhdy_j=h_n*n_y;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j);
           end %m
       end %i
   end %k
   end %l
end %e

%Boundary Conditions
for i=1:npoin
    j=iperiodic(i);
    if (i ~= j) 
        Mmatrix(i,:)=0;
        Mmatrix(i,i)=1;
        Lmatrix(i,:)=0;
        Lmatrix(i,i)=1;
    end
end


      
