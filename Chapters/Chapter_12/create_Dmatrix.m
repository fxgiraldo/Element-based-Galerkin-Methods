%---------------------------------------------------------------------%
%This function computes the 2D Differentiation Matrix on Quadrilaterals
%for DG in Weak Form
%Written by F.X. Giraldo on 8/24/2015
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Dmatrix(intma,jac,ksi_x,ksi_y,eta_x,eta_y,...
                                  psi,dpsi,npoin,nelem,ngl,nq,q)

%Initialize
rhs=zeros(npoin,2);
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
       
       %Interpolate onto Quadrature Points
       q_k=0; 
       for n=1:ngl
       for m=1:ngl
           jp=inode(m,n);
           q_k=q_k + psi(m,k)*psi(n,l)*q(jp);
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
           rhs(ip,1)=rhs(ip,1) - wq*dhdx_i*q_k;
           rhs(ip,2)=rhs(ip,2) - wq*dhdy_i*q_k;
       end %i
       end%j
       
   end %k
   end %l
end %e



      
