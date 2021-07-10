%---------------------------------------------------------------------%
%This function computes the RHS Volume contribution for the 1D SWE.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_volume(qp,qb,intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,gravity,delta_nl)

%Initialize
rhs=zeros(npoin,2);
hs_e=zeros(ngl,1);
U_e=zeros(ngl,1);
hb_e=zeros(ngl,1);
x=zeros(ngl,1);

%Integrate Divergence of Flux (Weak Form)
for e=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      ip=intma(i,e);
      hs_e(i)=qp(ip,1);
      U_e(i) =qp(ip,2);
      hb_e(i)=qb(ip);
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
         ip=intma(i,e);
         h_i=psi(i,l);
         dhdx_i=dpsi(i,l)*ksi_x;
         rhs(ip,1)=rhs(ip,1) + wq*dhdx_i*flux_h;
         rhs(ip,2)=rhs(ip,2) + wq*(dhdx_i*flux_U + h_i*gravity*hs_k*dhbdx_k);
      end %i
   end %l
end %e