%---------------------------------------------------------------------%
%This function computes the 2D Laplacian Matrix on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%Optimized by F.X. Giraldo on May 26, 2024 to reduce cost from 
% O(N^3d) to ~ O(N^2d) = N^d [ 2 N^2 ]= 2 N^{2d} for d=2
%In General we get O( d*N^{d + 2} )
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
   for j=1:nq
   for i=1:nq
       wq=wnq(i)*wnq(j)*jac(i,j,e);
       e_x=ksi_x(i,j,e);
       e_y=ksi_y(i,j,e);
       n_x=eta_x(i,j,e);
       n_y=eta_y(i,j,e);

       %Mass Matrix
       I=inode(i,j);
       Mmatrix(I,I)=Mmatrix(I,I) + wq;

       %Laplacian Matrix
       %Loop through I points = Rows of Matrix
       for ii=1:ngl

           %XI derivatives
           I=inode(ii,j); %j=ji
           h_e=dpsi(ii,i)*psi(j,j); %j=ji
           dhdx_i=h_e*e_x;
           dhdy_i=h_e*e_y;

           %Loop through J points = Columns of Matrix
           for ij=1:ngl
               J=inode(ij,j); %j=jj
               h_e=dpsi(ij,i)*psi(j,j); %j=jj
               dhdx_j=h_e*e_x;
               dhdy_j=h_e*e_y;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j);
               J=inode(i,ij); %i=ij and swap jj-> ij
               h_n=psi(i,i)*dpsi(ij,j);
               dhdx_j=h_n*n_x;
               dhdy_j=h_n*n_y;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j);
           end %ij

           %ETA derivatives
           I=inode(i,ii); %ii=i and swapped ij->ii
           h_n=psi(i,i)*dpsi(ii,j); %ii=i and swapped ij->ii
           dhdx_i=h_n*n_x;
           dhdy_i=h_n*n_y;
             
           for ij=1:ngl
               J=inode(ij,j); %jj=j
               h_e=dpsi(ij,i)*psi(j,j); %jj=j
               dhdx_j=h_e*e_x;
               dhdy_j=h_e*e_y;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j);
               J=inode(i,ij); %ij=i and swap jj-> ij
               h_n=psi(i,i)*dpsi(ij,j);
               dhdx_j=h_n*n_x;
               dhdy_j=h_n*n_y;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j);
           end %ij
       end %ii
   end %i
   end %j
end %e