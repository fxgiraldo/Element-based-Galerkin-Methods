%---------------------------------------------------------------------%
%This function computes the 2D Laplacian Matrix on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%Cost is ~ O(N^2d) = N^d[ 9 N^{d-1} ]= 9 N^5 for d=3
%In General we get O( d^d*N^{d+2} )
%---------------------------------------------------------------------%
function [Mmatrix,Lmatrix] = create_Global_Matrices_Inexact(intma,jac,wnq,ksi_x,ksi_y,ksi_z, ...
                             eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,psi,dpsi, ...
                             npoin,nelem,ngl,nq)

%Initialize
Mmatrix=zeros(npoin,npoin);
Lmatrix=zeros(npoin,npoin);
inode=zeros(ngl,ngl);

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
   for k=1:nq
   for j=1:nq
   for i=1:nq
       wq=wnq(i)*wnq(j)*wnq(k)*jac(i,j,k,e);
       e_x=ksi_x(i,j,k,e);
       e_y=ksi_y(i,j,k,e);
       e_z=ksi_z(i,j,k,e);
       n_x=eta_x(i,j,k,e);
       n_y=eta_y(i,j,k,e);
       n_z=eta_z(i,j,k,e);
       c_x=zeta_x(i,j,k,e);
       c_y=zeta_y(i,j,k,e);
       c_z=zeta_z(i,j,k,e);

       %Mass Matrix
       I=inode(i,j,k);
       Mmatrix(I,I)=Mmatrix(I,I) + wq;

       %Laplacian Matrix
       %Loop through I points = Rows of Matrix
       
       %XI derivatives
       for ii=1:ngl
      
           I=inode(ii,j,k); %j=ji, k=ki
           h_e=dpsi(ii,i)*psi(j,j)*psi(k,k); %j=ji, k=ki
           dhdx_i=h_e*e_x;
           dhdy_i=h_e*e_y;
           dhdz_i=h_e*e_z;
             
           %Loop through J points = Columns of Matrix
           for ij=1:ngl
               J=inode(ij,j,k); %j=jj, k=kj
               h_e=dpsi(ij,i)*psi(j,j)*psi(k,k);
               dhdx_j=h_e*e_x;
               dhdy_j=h_e*e_y;
               dhdz_j=h_e*e_z;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j + dhdz_i*dhdz_j);
           end %ij
       end %ii

       %ETA derivatives
       for ji=1:ngl
      
           I=inode(i,ji,k); %i=ii, k=ki
           h_n=psi(i,i)*dpsi(ji,j)*psi(k,k); %i=ii, k=ki
           dhdx_i=h_n*n_x;
           dhdy_i=h_n*n_y;
           dhdz_i=h_n*n_z;
             
           %Loop through J points = Columns of Matrix
           for jj=1:ngl
               J=inode(i,jj,k); %i=ij, k=kj
               h_n=psi(i,i)*dpsi(jj,j)*psi(k,k);
               dhdx_j=h_n*n_x;
               dhdy_j=h_n*n_y;
               dhdz_j=h_n*n_z;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j + dhdz_i*dhdz_j);
           end %jj
       end %ji

       %ZETA derivatives
       for ki=1:ngl
      
           I=inode(i,j,ki); %i=ii, j=ji
           h_c=psi(i,i)*psi(j,j)*dpsi(ki,k); %i=ii, j=ji
           dhdx_i=h_c*c_x;
           dhdy_i=h_c*c_y;
           dhdz_i=h_c*c_z;
             
           %Loop through J points = Columns of Matrix
           for kj=1:ngl
               J=inode(i,j,kj); %i=ij, j=jj
               h_c=psi(i,i)*psi(j,j)*dpsi(kj,k);
               dhdx_j=h_c*c_x;
               dhdy_j=h_c*c_y;
               dhdz_j=h_c*c_z;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j + dhdz_i*dhdz_j);
           end %kj
       end %ki

   end %i
   end %j
   end %k
end %e

      
