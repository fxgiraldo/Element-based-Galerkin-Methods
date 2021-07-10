%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Naval Research Laboratory 
%           Monterey, CA 93943-5502
%---------------------------------------------------------------------%
function rhs = create_rhs_volume(qp,intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,alpha)

%Initialize
rhs=zeros(npoin,2);
U_e=zeros(ngl,1);
x=zeros(ngl,1);

%Integrate Divergence of Flux (Weak Form)
for e=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      ip=intma(i,e);
      U_e(i)=qp(ip,1);
      x(i)  =coord(i,e);
   end
   
   %Jacobians
   dx=x(ngl)-x(1);
   jac=dx/2;
   ksi_x=2/dx;
   
   %LGL Integration
   for k=1:nq
      wq=wnq(k)*jac;
            
      %Form Derivative
      U_k=0;
      f_k=0;
      dUdx_k=0;
      for j=1:ngl
         h_j=psi(j,k);
         dhdx_j=dpsi(j,k)*ksi_x;
         U_k=U_k + h_j*U_e(j);
         f_k=f_k + h_j*(0.5*U_e(j)^2);
         dUdx_k=dUdx_k + dhdx_j*(U_e(j));
      end
      
      %Form RHS
      for i=1:ngl
          ip=intma(i,e);
          h_i=psi(i,k);
          U_i=U_e(i);
          dhdx_i=dpsi(i,k)*ksi_x;
          dhUdx_k=dhdx_i*U_k;
          %Both RHS below give the same answer = equivalent for Burgers
          rhs(ip,1)=rhs(ip,1) + wq*( dhdx_i*f_k ...
                              - (1.0-alpha)*h_i*(0.5*U_k)*dUdx_k );
%           rhs(ip,1)=rhs(ip,1) + wq*( alpha*dhdx_i*f_k ...
%                               + (1.0-alpha)*(dhUdx_k*(0.5*U_k)-h_i*(0.5*U_k)*dUdx_k) );
      end %i
   end %l
end %e