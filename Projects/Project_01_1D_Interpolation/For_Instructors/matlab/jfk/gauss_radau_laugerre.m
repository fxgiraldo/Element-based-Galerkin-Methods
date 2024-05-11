function [xgr,wgr] = gauss_radau_laugerre(P,scale)
% [xgl,wgl] = gauss_radau_laugerre(P,scale)
% Computes the the Gauss Radau Laguerre qudrature points and weights
% for the Laguerre function L_i (x) if scale = false.
% If scale =true, then returns points/weights
% for the scaled Laguerre function exp(-x/2) * L_i (x), where 
% L_i(x) is the i-th Laguerre polynomial

% See formulas in Valenciano and Chaplain, Mathematical and Computer 
% Modelling 41 (2005) 1171-1192

% James F. Kelly
% 3 February 2023

% Use matlab root finder as a first guess
Pp1 = P+1;
xgr =  roots(polyder(LaguerrePoly(Pp1)));  
xgr = flipud([xgr; 0]);
xgr = real(xgr);
    
Lk = LaguerrePoly(P);
Lkx = polyval(Lk,xgr);
wgr = 1./(Pp1.*Lkx.^2);
if(scale)
   wgr = exp(xgr).*wgr; 
end
  
end
