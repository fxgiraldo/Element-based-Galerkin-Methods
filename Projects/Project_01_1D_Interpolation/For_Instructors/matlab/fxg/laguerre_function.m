%---------------------------------------------------------------------%
%This code computes the Laguerre Functions and its 1st derivatives
%Written by F.X. Giraldo on 11/20/2023
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [L0] = laguerre_function(p,x)

L0=exp(-0.5*x);
L1=exp(-0.5*x)*(1-x);

for i=2:p
    L2=L1; L1=L0;
    L0=((2*i-1-x)/i)*L1 - ((i-1)/i)*L2;
end
