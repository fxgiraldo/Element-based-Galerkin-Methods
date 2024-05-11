function [psi,dpsi] = lagrange_basis_laguerre_functions(xgr)
%[psi,dpsi] = lagrange_basis_laguerre_functions(xgl)
% See Eq. (3.17) in Shen 2000

nbasis = length(xgr);
N = nbasis -1;
Np1 = N+1;

psi = eye(nbasis);
dpsi = zeros(nbasis);
dpsi(1,1)=-(N +1)./2.0;

Lkxi = scaled_laguerre(xgr,Np1);


for i = 1:nbasis
    xi = xgr(i);
    for j = 1:nbasis
        xj = xgr(j);
        if(i ~= j)
           %dpsi(j,i) = scaled_laguerre(xi,Np1)/(scaled_laguerre(xj,Np1).*(xi -xj));
           dpsi(j,i) = Lkxi(i)/(Lkxi(j).*(xi-xj));
        end
    end
end

end


function y = scaled_laguerre(x,n)
Lkx = real(laguerreL(n,x));
y = exp(-x/2).*Lkx;
end