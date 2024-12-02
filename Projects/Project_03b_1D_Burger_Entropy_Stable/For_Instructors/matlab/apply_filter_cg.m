%---------------------------------------------------------------------%
%This function applies the Legendre Filter.
%Written by F.X. Giraldo on 7/2007
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_filter_cg(qp,fq,intma,coord,npoin,nelem,ngl,nq,wnq,psi,iperiodic)

%Initialize
rhs=zeros(npoin,2);
q_e=zeros(ngl,1);
inode=zeros(ngl,1);
x=zeros(ngl,1);

%Integrate Divergence of Flux
for e=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      I=intma(i,e);
      inode(i)=iperiodic(I);
      x(i)=coord(i,e);
      q_e(i)=qp(I,1);
   end

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
         I=inode(i);
         h_i=psi(i,l);
         rhs(I,1)=rhs(I,1) + wq*h_i*qf_k;
      end %i
   end %l
end %ie

%Periodicity
% if (iperiodic(npoin) ~= npoin)
%     rhs(npoin)=0;
% end
