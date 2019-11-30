%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs(qp,coord,nelem,ngl,nq,wnq,psi,dpsi,u)

%Initialize
rhs=zeros(ngl,nelem);
q_e=zeros(ngl,1);
x=zeros(ngl,1);

%Integrate Divergence of Flux
for e=1:nelem
   
   %Store Coordinates
   q_e(:)=qp(:,e);
   x(:)=coord(:,e);
   
   %Jacobians
   dx=x(ngl)-x(1);
   jac=dx/2;
   ksi_x=2/dx;
   
   %LGL Integration
   for l=1:nq
      wq=wnq(l)*jac;
            
      %Form Derivative
      q_k=0;
      q_ksi=0;
      for j=1:ngl
          q_k=q_k + psi(j,l)*q_e(j)*u;
          q_ksi=q_ksi + dpsi(j,l)*q_e(j);
      end
      q_x=q_ksi*ksi_x;
      
      %Form RHS
      for i=1:ngl
         dhdx_i=dpsi(i,l)*ksi_x;
%          rhs(i,e)=rhs(i,e) + wq*dhdx_i*q_k - 0*wq*visc*dhdx_i*q_x;
         rhs(i,e)=rhs(i,e) + wq*dhdx_i*q_k;

      end %i
      
   end %l
end %e
