%---------------------------------------------------------------------%
%This function computes the 2D Laplacian Matrix on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Lmatrix = create_Lmatrix(intma,jac,wnq,ksi_x,ksi_y,ksi_z,...
                   eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,psi,dpsi,...
                   npoin,nelem,ngl,nq)

%Initialize
Lmatrix=zeros(npoin,npoin);
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
       e_x=ksi_x(l,m,n,e);
       e_y=ksi_y(l,m,n,e);
       e_z=ksi_z(l,m,n,e);
       n_x=eta_x(l,m,n,e);
       n_y=eta_y(l,m,n,e);
       n_z=eta_z(l,m,n,e);
       c_x=zeta_x(l,m,n,e);
       c_y=zeta_y(l,m,n,e);
       c_z=zeta_z(l,m,n,e);
       
       %Loop through I points
       for k=1:ngl
       for j=1:ngl
       for i=1:ngl
           I=inode(i,j,k);
           h_e=dpsi(i,l)*psi(j,m)*psi(k,n);
           h_n=psi(i,l)*dpsi(j,m)*psi(k,n);
           h_c=psi(i,l)*psi(j,m)*dpsi(k,n);
           dhdx_i=h_e*e_x + h_n*n_x + h_c*c_x;
           dhdy_i=h_e*e_y + h_n*n_y + h_c*c_y;
           dhdz_i=h_e*e_z + h_n*n_z + h_c*c_z;           
           
           %Loop through J points
           for kj=1:ngl
           for jj=1:ngl
           for ij=1:ngl
               J=inode(ij,jj,kj);
               h_e=dpsi(ij,l)*psi(jj,m)*psi(kj,n);
               h_n=psi(ij,l)*dpsi(jj,m)*psi(kj,n);
               h_c=psi(ij,l)*psi(jj,m)*dpsi(kj,n);
               dhdx_j=h_e*e_x + h_n*n_x + h_c*c_x;
               dhdy_j=h_e*e_y + h_n*n_y + h_c*c_y;
               dhdz_j=h_e*e_z + h_n*n_z + h_c*c_z;
               Lmatrix(I,J)=Lmatrix(I,J) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j + dhdz_i*dhdz_j);
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