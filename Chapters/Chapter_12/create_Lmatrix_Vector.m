%---------------------------------------------------------------------%
%This function computes the 2D Laplacian Matrix on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Rvector = create_Lmatrix_Vector(intma,jac,ksi_x,ksi_y,eta_x,...         
                   eta_y,psi,dpsi,npoin,nelem,ngl,nq,q)

%Initialize
Rvector=zeros(npoin,1);
inode=zeros(ngl,ngl);

for e=1:nelem
   
   %Store Coordinates
   for j=1:ngl
   for i=1:ngl
      ip=intma(e,i,j);
      inode(i,j)=ip;
   end %i
   end %j
   
   %Do LGL Integration
   for l=1:nq
   for k=1:nq
       wq=jac(e,k,l);
       e_x=ksi_x(e,k,l);
       e_y=ksi_y(e,k,l);
       n_x=eta_x(e,k,l);
       n_y=eta_y(e,k,l);
       
       %Interpolate Derivatives onto Quadrature Points
       q_x=0; q_y=0;
       for n=1:ngl
       for m=1:ngl
           jp=inode(m,n);
           h_e=dpsi(m,k)*psi(n,l);
           h_n=psi(m,k)*dpsi(n,l);
           dhdx_j=h_e*e_x + h_n*n_x;
           dhdy_j=h_e*e_y + h_n*n_y;
           q_x=q_x + dhdx_j*q(jp);
           q_y=q_y + dhdy_j*q(jp);
       end %m
       end %n
           
       %Loop through I points
       for j=1:ngl
       for i=1:ngl
           ip=inode(i,j);
           h_e=dpsi(i,k)*psi(j,l);
           h_n=psi(i,k)*dpsi(j,l);
           dhdx_i=h_e*e_x + h_n*n_x;
           dhdy_i=h_e*e_y + h_n*n_y;
           Rvector(ip)=Rvector(ip) - wq*(dhdx_i*q_x + dhdy_i*q_y);
       end %i
       end%j
   end %k
   end %l
end %e

      
