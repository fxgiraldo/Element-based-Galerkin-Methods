%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Lmatrix_IBP_General(intma,jac,npoin,nelem,ngl,wnq,dpsi,q,visc_elem)

%Initialize
rhs=zeros(npoin,2);

for e=1:nelem
   
    h_e=zeros(ngl(e),1);
    U_e=zeros(ngl(e),1);
    visc_h=zeros(ngl(e),1);
    visc_u=zeros(ngl(e),1);

   %Store Coordinates
   for i=1:ngl(e)
      I=intma(i,e);
      h_e(i)=q(I,1);
      U_e(i)=q(I,2);
      visc_h(i)=visc_elem(I,1);
      visc_u(i)=visc_elem(I,2);
   end
   
   %Jacobians
   jac_e=jac(e);
   ksi_x=1.0/jac_e;
   
   for i=1:ngl(e)
      wq=wnq(i,e)*jac_e;
                  
      %Form Derivative
      h_x=0;
      U_x=0;
      for j=1:ngl(e)
          h_x=h_x + dpsi(j,i,e)*ksi_x*h_e(j);
          U_x=U_x + dpsi(j,i,e)*ksi_x*U_e(j);
      end
      
      %Form RHS
      for j=1:ngl(e)
         J=intma(j,e);
         dhdx_ji=dpsi(j,i,e)*ksi_x;
         rhs(J,1)=rhs(J,1) - wq*visc_h(i)*dhdx_ji*h_x;
         rhs(J,2)=rhs(J,2) - wq*visc_u(i)*dhdx_ji*U_x;
      end %i
   end %l
   
end %e

