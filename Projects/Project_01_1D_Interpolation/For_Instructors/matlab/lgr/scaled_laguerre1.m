function y = scaled_laguerre1(n,x)
% Scaled Laguerre functions that satsify a zero Dirichlet BC at left

Lk = laguerreL(n,x);
%if (n < 2) 
    y = exp(-x/2).*Lk;
%else
%    Lk2 = laguerreL(n-2,x);
%    y = exp(-x/2).*(Lk-Lk2);
%end

end