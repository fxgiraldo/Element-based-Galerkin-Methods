%---------------------------------------------------------------------%
%This code computes the Laguerre Polynomials and its 1st derivatives
%Written by F.X. Giraldo on 11/20/2023
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [L0,L0_1] = laguerre_poly(p,x)

L0=1;
L1=1-x;
L0_1=0;
L1_1=-1;

for i=2:p
    L2=L1;L2_1=L1_1;
    L1=L0;L1_1=L0_1;
    L0=(2*i-1-x)/i*L1 - (i-1)/i*L2;
    L0_1=(2*i-x)/i*L1_1 - L2_1;
end
