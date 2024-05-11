%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Naval Research Laboratory 
%           Monterey, CA 93943-5502
%---------------------------------------------------------------------%
function rhs = compute_volume_integral(qp,intma,coord,npoin,nelem,ngl,nq,wnq,psi)

%Initialize
rhs=zeros(npoin,2);
h_e=zeros(ngl,1);
U_e=zeros(ngl,1);
x=zeros(ngl,1);

%Integrate Divergence of Flux (Weak Form)
for e=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      ip=intma(i,e);
      h_e(i)=qp(ip,1);
      U_e(i)=qp(ip,2);
      x(i)=coord(i,e);
   end
   
   %Jacobians
   dx=x(ngl)-x(1);
   jac=dx/2;
   
   %LGL Integration
   for l=1:nq
      wq=wnq(l)*jac;
                  
      %Interpolate onto Quadrature Points
      h_k=0;
      U_k=0;
      for j=1:ngl
         h_j=psi(j,l);
         h_k=h_k + h_e(j)*h_j;
         U_k=U_k + U_e(j)*h_j;
      end

      %Form RHS
      for i=1:ngl
         ip=intma(i,e);
         h_i=psi(i,l);
         rhs(ip,1)=rhs(ip,1) + wq*h_i*h_k;
         rhs(ip,2)=rhs(ip,2) + wq*h_i*U_k;
      end %i
   end %l
end %e