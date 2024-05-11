function [psi,dpsi] = lagrange_basis_laguerre_v2(P,Q,xgl,xnq)
%[psi,dpsi] = lagrange_basis_laguerre_functions(xgl)
% See Eq. (3.17) in Shen 2000

% nbasis = length(xgr);
% N = nbasis -1;
% Np1 = N+1;

% psi = eye(nbasis);
% dpsi = zeros(nbasis);
% dpsi(1,1)=-(N +1)./2.0;

psi=zeros(P,Q);
dpsi=zeros(P,Q);

%Build Basis Function
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
      end
   end
end

%Build Derivative of Scaled Basis Function
Lkxi = scaled_laguerre(xgl,P);
for i = 1:Q
    xi = xnq(i);
    for j = 1:P
        xj = xgl(j);
        if(i ~= j)
           dpsi(j,i) = Lkxi(i)/(Lkxi(j).*(xi-xj));
        end
    end
end

end


function y = scaled_laguerre(x,n)
Lkx = real(laguerreL(n,x));
y = exp(-x/2).*Lkx;
end