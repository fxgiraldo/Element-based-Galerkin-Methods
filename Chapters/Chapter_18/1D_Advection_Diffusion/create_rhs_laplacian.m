%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Naval Postgraduate School 
%           Monterey, CA 93943-5502
%---------------------------------------------------------------------%
function rhs = create_rhs_laplacian(rhs,qp,coord,nelem,ngl,nq,wnq,psi,dpsi)

%Initialize
%rhs=zeros(ngl,nelem);
q_e=zeros(ngl,1);
x=zeros(ngl,1);

%Integrate Divergence of Flux
for ie=1:nelem
   
   %Store Coordinates
   q_e(:)=qp(:,ie);
   x(:)=coord(:,ie);
   
   %Jacobians
   dx=x(ngl)-x(1);
   jac=dx/2;
   ksi_x=2/dx;
   
   %LGL Integration
   for l=1:nq
      wq=wnq(l)*jac;
            
      %Form Derivative
      q_k=0;
      for j=1:ngl
          q_k=q_k + psi(j,l)*q_e(j)*u;
      end
      
      %Form RHS
      for i=1:ngl
         dhdx_i=dpsi(i,l)*ksi_x;
         rhs(i,ie)=rhs(i,ie) + wq*dhdx_i*q_k;
      end %i
   end %l
end %ie
