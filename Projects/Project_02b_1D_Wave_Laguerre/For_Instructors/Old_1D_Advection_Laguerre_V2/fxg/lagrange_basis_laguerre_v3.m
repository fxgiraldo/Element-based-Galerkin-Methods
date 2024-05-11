%---------------------------------------------------------------------%
%This code computes the cardinal functions based on Laguerre polynomials  
% at the sample points XNQ. It uses the standard Lagrange polynomial procedure 
% and then adjusts to form the cardinal functions based on the Laguerre
% functions.
%Written by F.X. Giraldo on 11/2023
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [psiL,dpsiL,xnq,wnq] = lagrange_basis_laguerre_v3(P,Q,xgl)

[xnq,wnq]= laguerre_gauss_radau_eigenvalue(Q,true); %true=scale weights by exp(x)

%Lagrange Polynomials based on roots of Orthogonal Polynomails (Legendre or Laguerre)
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
               if (k ~=i && k ~=j)
                  ddpsi=ddpsi*(xl-xk)/(xi-xk);
               end
            end
            dpsi(i,l)=dpsi(i,l) + ddpsi/(xi-xj);
         end
      end
   end
end

%Using Laguerre Polynomial cardinal functions (psi) and its derivative
%(dpsi), build for corresponding Laguerre Functions.
for i=1:P
    for j=1:Q
        factor=exp(-0.5*(xnq(j)-xgl(i)));
        psiL(i,j)=factor*psi(i,j);
        dpsiL(i,j)=factor*( dpsi(i,j) - 0.5*psi(i,j) );
    end
end
