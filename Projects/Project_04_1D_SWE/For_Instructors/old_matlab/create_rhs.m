%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Naval Research Laboratory 
%           Monterey, CA 93943-5502
%---------------------------------------------------------------------%
function rhs = create_rhs_volume(qp,qb,coord,nelem,ngl,nq,wnq,psi,dpsi,gravity,delta_nl)

%Initialize
rhs=zeros(2,ngl,nelem);
hs_e=zeros(ngl,1);
U_e=zeros(ngl,1);
hb_e=zeros(ngl,1);
x=zeros(ngl,1);

%Integrate Divergence of Flux (Weak Form)
for e=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      hs_e(i)=qp(1,i,e);
      U_e(i) =qp(2,i,e);
      hb_e(i)=qb(i,e);
      x(i)   =coord(i,e);
   end
   
   %Jacobians
   dx=x(ngl)-x(1);
   jac=dx/2;
   ksi_x=2/dx;
   
   %LGL Integration
   for l=1:nq
      wq=wnq(l)*jac;
            
      %Form Derivative
      hs_k=0;
      U_k=0;
      hb_k=0;
      dhbdx_k=0;
      for j=1:ngl
         h_j=psi(j,l);
         dhdx_j=dpsi(j,l)*ksi_x;
         U_k=U_k + U_e(j)*h_j;
         hb_k =hb_k  + hb_e(j)*h_j;
         hs_k=hs_k + hs_e(j)*h_j;
         dhbdx_k=dhbdx_k + hb_e(j)*dhdx_j;
      end
      
      h_k=hs_k + hb_k;
      flux_h=U_k;
      flux_U=( U_k^2/(h_k) + 0.5*gravity*hs_k^2 )*delta_nl + gravity*hb_k*hs_k;
      
      %Form RHS
      for i=1:ngl
         h_i=psi(i,l);
         dhdx_i=dpsi(i,l)*ksi_x;
         rhs(1,i,e)=rhs(1,i,e) + wq*dhdx_i*flux_h;
         rhs(2,i,e)=rhs(2,i,e) + wq*(dhdx_i*flux_U + h_i*gravity*hs_k*dhbdx_k);
      end %i
   end %l
end %e