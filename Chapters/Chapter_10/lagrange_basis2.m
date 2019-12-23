%---------------------------------------------------------------------%
%This code computes the Legendre cardinal functions at the Quadrature 
%points
%Written by F.X. Giraldo on 4/2000
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [psi,dpsi] = lagrange_basis2(P,xgl,Ns,xs)

psi=zeros(P,Ns);

for l=1:Ns
    xl=xs(l);
   
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
         end %if
         ddpsi=1;
         if (j ~= i)
            for k=1:P
               xk=xgl(k);
               %Derivative of Basis
               if (k ~=i & k ~=j)
                  ddpsi=ddpsi*(xl-xk)/(xi-xk);
               end %if
            end %k
            dpsi(i,l)=dpsi(i,l) + ddpsi/(xi-xj);
         end %if
      end %j
   end %i
end %l


      
