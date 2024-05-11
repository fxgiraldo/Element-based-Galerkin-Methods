function [xgr,wgr] = gauss_radau_laguerre_eig(P,scale)
% Uses the eigenvalue method to get a first guess for the LGR nodes, and 
% then refines via Newton's method
% See Chapter 3 in thesis by Tianyao Yue: "Spectral Element Method for
% Pricing European Options"
% Note: Use the internal 

% James F. Kelly
% 24 February 2023


Pp1 = P +1;
n = 0:P;

%See matrix 3.243 and 3.244 in Yue's thesis
an = 2.*n + 1;
an(P+1) = P;
bn = 1:P;
J = diag(an) + diag(bn,1) + diag(bn,-1);
xi = eig(J);
xgr = xi;           %Initial guess

ngr = length(xgr);
thresh = 1e-10;

 for k = 1:ngr
       x0 = xgr(k);
       diff1 = 1.0;
      while(diff1 > thresh)
          x1 = x0 + (laguerreL(Pp1,x0) - laguerreL(P,x0))/laguerreL(P,x0);
          diff1 = abs(x1 -x0);
          x0 = x1;
      end
      xgr(k) = x1;
 end
 xgr(1) = 0;
 
 
Lkx = laguerreL(P,xgr);
wgr = 1./(Pp1.*Lkx.^2);
if(scale)
   wgr = exp(xgr).*wgr; 
end

end
