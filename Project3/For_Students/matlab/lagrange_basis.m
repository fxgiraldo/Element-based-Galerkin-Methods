%---------------------------------------------------------------------%
%This code computes the Legendre cardinal functions at the Quadrature 
%points
%Written by F.X. Giraldo on 4/2000
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [psi,dpsi] = lagrange_basis(P,Q,xgl,xnq)

%Get Quadrature roots
%[xnq,wnq] = legendre_gauss_lobatto(Q);

%Perform Quadrature
for l=1:Q
   xl=xnq(l);
   
   %Construct Basis
   for i=1:P
      xi=xgl(i);
      psi(i,l)=1;
      dpsi(i,l)=0;
      for j=1:P
         xj=xgl(j);
         %Basis
         if (j ~= i)
            psi(i,l)=psi(i,l)*(xl-xj)/(xi-xj);
         end
         ddpsi=1;
         if (j ~= i)
            for k=1:P
               xk=xgl(k);
               %Derivative of Basis
               if (k ~=i & k ~=j)
                  ddpsi=ddpsi*(xl-xk)/(xi-xk);
               end
            end
            dpsi(i,l)=dpsi(i,l) + ddpsi/(xi-xj);
         end
      end
   end
end
