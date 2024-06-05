%---------------------------------------------------------------------%
%This function computes the 2D Laplacian Matrix on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%Cost is O( N^{3d} )
%---------------------------------------------------------------------%
function [Mmatrix,Lmatrix] = create_Global_Matrices_Exact(intma,jac,wnq,ksi_x,ksi_y,ksi_z, ...
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
       
       %Loop through I points
       for ki=1:ngl
       for ji=1:ngl
       for ii=1:ngl
           I=inode(ii,ji,ki);
           h_i=psi(ii,i)*psi(ji,j)*psi(ki,k);
           h_e=dpsi(ii,i)*psi(ji,j)*psi(ki,k);
           h_n=psi(ii,i)*dpsi(ji,j)*psi(ki,k);
           h_c=psi(ii,i)*psi(ji,j)*dpsi(ki,k);
           dhdx_i=h_e*e_x + h_n*n_x + h_c*c_x;
           dhdy_i=h_e*e_y + h_n*n_y + h_c*c_y;
           dhdz_i=h_e*e_z + h_n*n_z + h_c*c_z;
             
           %Loop through J points
           for kj=1:ngl
           for jj=1:ngl
           for ij=1:ngl
               J=inode(ij,jj,kj);
               h_j=psi(ij,i)*psi(jj,j)*psi(kj,k);
               h_e=dpsi(ij,i)*psi(jj,j)*psi(kj,k);
               h_n=psi(ij,i)*dpsi(jj,j)*psi(kj,k);
               h_c=psi(ij,i)*psi(jj,j)*dpsi(kj,k);
               dhdx_j=h_e*e_x + h_n*n_x + h_c*c_x;
               dhdy_j=h_e*e_y + h_n*n_y + h_c*c_y;
               dhdz_j=h_e*e_z + h_n*n_z + h_c*c_z;
               Mmatrix(I,J)=Mmatrix(I,J) + wq*h_i*h_j;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j + dhdz_i*dhdz_j);
           end %ij
           end %jj
           end %kj
       end %ii
       end %ji
       end %ki
   end %i
   end %j
   end %k
end %e

      
