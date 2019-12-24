%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Maths
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs(qp,qb,coord,nelem,ngl,nq,wnq,psi,dpsi,h_eps,gravity,delta_nl)

%Initialize
rhs=zeros(2,ngl,nelem);
ps_e=zeros(ngl,1);
pu_e=zeros(ngl,1);
u_e=zeros(ngl,1);
b_e=zeros(ngl,1);
x=zeros(ngl,1);

%Integrate Divergence of Flux (Weak Form)
for ie=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      ps_e(i)=qp(1,i,ie);
      pu_e(i)=qp(2,i,ie);
      b_e(i) =qb(i,ie);
      x(i)   =coord(i,ie);
   end
   
   %Jacobians
   dx=x(ngl)-x(1);
   jac=dx/2;
   ksi_x=2/dx;
   
   %LGL Integration
   for l=1:nq
      wq=wnq(l)*jac;
            
      %Form Derivative
      ps_k=0;
      pu_k=0;
      b_k=0;
      dbdx=0;
      for j=1:ngl
         h_k=psi(j,l);
         dhdx=dpsi(j,l)*ksi_x;
         pu_k=pu_k + pu_e(j)*h_k;
         b_k =b_k  + b_e(j)*h_k;
         ps_k=ps_k + ps_e(j)*h_k;
         dbdx=dbdx + b_e(j)*dhdx;
      end
      
      %-----New
      p_k=ps_k + b_k;
      if (p_k <= h_eps)
          p_k=h_eps;
          %ps_k=p_k - b_k;
          pu_k=0;
      end
      %-----New
      
      flux_p=pu_k;
      flux_pu=( pu_k^2/(p_k) + 0.5*gravity*ps_k^2 )*delta_nl + gravity*b_k*ps_k;
   
      %Form RHS
      for i=1:ngl
         h_i=psi(i,l);
         dhdx_i=dpsi(i,l)*ksi_x;
         rhs(1,i,ie)=rhs(1,i,ie) + wq*dhdx_i*flux_p;
         rhs(2,i,ie)=rhs(2,i,ie) + wq*(dhdx_i*flux_pu + h_i*gravity*ps_k*dbdx);

      end %i
   end %l
end %ie
