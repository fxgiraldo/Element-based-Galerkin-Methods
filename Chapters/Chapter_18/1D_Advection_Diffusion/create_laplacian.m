%---------------------------------------------------------------------%
%This function computes the Laplacian in 1D.
%Written by F.X. Giraldo on 9/2012
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_laplacian(qp,coord,nelem,ngl,nq,wnq,dpsi)

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
      q_ksi=0;
      for j=1:ngl
          q_ksi=q_ksi + dpsi(j,l)*q_e(j);
      end
      q_x=q_ksi*ksi_x;
      
      %Form RHS
      for i=1:ngl
         dhdx_i=dpsi(i,l)*ksi_x;
         rhs(i,e)=rhs(i,e) - wq*dhdx_i*q_x;
      end %i
   end %l
end %e
