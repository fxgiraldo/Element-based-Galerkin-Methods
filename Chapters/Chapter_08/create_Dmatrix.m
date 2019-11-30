%---------------------------------------------------------------------%
%This function computes the Volume Integral Part 
%of the Differentiation Matrix.
%Written by F.X. Giraldo on 8/17/2015
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Dmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,q)

%Initialize
rhs=zeros(npoin,1);
x=zeros(ngl,1);

for e=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      x(i)=coord(i,e);
   end
   
   dx=x(ngl)-x(1);
   jac=dx/2;
   ksi_x=2/dx;
   
   for l=1:nq
      wq=wnq(l)*jac;
            
      %interpolate q
      q_e=0;
      for j=1:ngl
          jp=intma(j,e);
          q_e=q_e + psi(j,l)*q(jp);
      end
      
      %Form RHS
      for i=1:ngl
          dhdx_i=dpsi(i,l)*ksi_x;
          ip=intma(i,e);
          rhs(ip)=rhs(ip) - wq*dhdx_i*q_e;
      end %i
   end %l
      
end %e