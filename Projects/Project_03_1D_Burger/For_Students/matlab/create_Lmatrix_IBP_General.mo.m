%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Lmatrix_IBP_General(intma,coord,npoin,nelem,ngl,nq,wnq,dpsi,q,visc_elem)

%Initialize
rhs=zeros(npoin,2);
x=zeros(ngl,1);
h_e=zeros(ngl,1);
U_e=zeros(ngl,1);
visc_h=zeros(ngl,1);
visc_u=zeros(ngl,1);

for e=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      ip=intma(i,e);
      h_e(i)=q(ip,1);
      U_e(i)=q(ip,2);
      x(i)=coord(i,e);
      visc_h(i)=visc_elem(ip,1);
      visc_u(i)=visc_elem(ip,2);
   end
   
   dx=x(ngl)-x(1);
   jac=dx/2;
   ksi_x=2/dx;
   
   for l=1:nq
      wq=wnq(l)*jac;
            
      %Form Derivative
      h_x=0;
      U_x=0;
      for j=1:ngl
          h_x=h_x + dpsi(j,l)*ksi_x*h_e(j);
          U_x=U_x + dpsi(j,l)*ksi_x*U_e(j);
      end
      
      %Form RHS
      for i=1:ngl
         ip=intma(i,e);
         dhdx_i=dpsi(i,l)*ksi_x;
         rhs(ip,1)=rhs(ip,1) - wq*visc_h(i)*dhdx_i*h_x;
         rhs(ip,2)=rhs(ip,2) - wq*visc_u(i)*dhdx_i*U_x;
      end %i
   end %l
   
end %e

