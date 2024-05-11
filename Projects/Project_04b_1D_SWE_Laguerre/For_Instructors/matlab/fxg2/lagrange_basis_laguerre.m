%---------------------------------------------------------------------%
%This code computes the Legendre cardinal functions at the Quadrature 
%points
%Written by F.X. Giraldo on 4/2000
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [psiL,dpsiL] = lagrange_basis_laguerre(P,Q,xgl)

%Get Quadrature roots
xnq = xgl;

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
               if (k ~=i & k ~=j)
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
        factor=exp(-0.5*xnq(j))/exp(-0.5*xgl(i));
        psiL(i,j)=factor*psi(i,j);
        delta=0;
        if (i == j) 
            delta=1;
        end
        dpsiL(i,j)=factor*( dpsi(i,j) - 0.5*delta);
    end
end
