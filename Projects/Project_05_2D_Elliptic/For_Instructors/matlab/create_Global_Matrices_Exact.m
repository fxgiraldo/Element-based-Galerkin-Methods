%---------------------------------------------------------------------%
%This function computes the 2D Laplacian Matrix on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%Cost is O(N^3d)
%---------------------------------------------------------------------%
function [Mmatrix,Lmatrix] = create_Global_Matrices_Exact(intma,jac,wnq,ksi_x,ksi_y,eta_x,...         
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
       
       %Loop through I points = Rows of Matrix
       for ji=1:ngl
       for ii=1:ngl
           I=inode(ii,ji);
           h_i=psi(ii,i)*psi(ji,j);
           h_e=dpsi(ii,i)*psi(ji,j);
           h_n=psi(ii,i)*dpsi(ji,j);
           dhdx_i=h_e*e_x + h_n*n_x;
           dhdy_i=h_e*e_y + h_n*n_y;
             
           %Loop through J points = Cols of Matrix
           for jj=1:ngl
           for ij=1:ngl
               J=inode(ij,jj);
               h_j=psi(ij,i)*psi(jj,j);
               h_e=dpsi(ij,i)*psi(jj,j);
               h_n=psi(ij,i)*dpsi(jj,j);
               dhdx_j=h_e*e_x + h_n*n_x;
               dhdy_j=h_e*e_y + h_n*n_y;
               Mmatrix(I,J)=Mmatrix(I,J) + wq*h_i*h_j;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j);
           end %jj
           end %ij
       end %ii
       end%ji
   end %i
   end %j
end %e

      
