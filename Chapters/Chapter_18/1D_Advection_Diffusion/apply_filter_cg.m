%---------------------------------------------------------------------%
%This function applies the Legendre Filter.
%Written by F.X. Giraldo on 7/2007
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_filter_cg(qp,fq,coord,nelem,ngl,nq,wnq,psi)

%Initialize
rhs=zeros(ngl,nelem);
q_e=zeros(ngl,1);
qf=zeros(nq,1);
x=zeros(ngl,1);

%Integrate Divergence of Flux
for e=1:nelem
   
   %Store Coordinates
  x(:)=coord(:,e);
  q_e(:)=qp(:,e);
  
   %Form Filtered Variable
   qf=fq*q_e;

   %Jacobians
   dx=x(ngl)-x(1);
   jac=dx/2;
   ksi_x=2/dx;
   
   %LGL Integration
   for l=1:nq
      wq=wnq(l)*jac;
            
      %Interpolate onto Quadrature Points
      qf_k=qf(l);

      %Form RHS
      for i=1:ngl
         h_i=psi(i,l);
         rhs(i,e)=rhs(i,e) + wq*h_i*qf_k;
      end %i
   end %l
end %e
