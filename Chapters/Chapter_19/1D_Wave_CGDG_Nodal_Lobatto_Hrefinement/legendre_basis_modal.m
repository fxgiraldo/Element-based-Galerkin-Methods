%---------------------------------------------------------------------%
%This code computes the Legendre Polynomials and its 1st and 2nd derivatives
%Written by F.X. Giraldo on 4/2000
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [L,dL] = legendre_basis_modal(P,Q,xnq)

L=zeros(P,Q);
dL=zeros(P,Q);
p=P-1;

for k=1:Q
    x=xnq(k);
    
    L1=0;L1_1=0;L1_2=0;
    L0=1;L0_1=0;L0_2=0;
    L(1,k)=L0;
    dL(1,k)=L0_1;
    for i=1:p
       L2=L1;L2_1=L1_1;L2_2=L1_2;
       L1=L0;L1_1=L0_1;L1_2=L0_2;
       a=(2*i-1)/i;
       b=(i-1)/i;
       L0=a*x*L1 - b*L2;
       L0_1=a*(L1 + x*L1_1) - b*L2_1;
       L0_2=a*(2*L1_1 + x*L1_2) - b*L2_2;
       L(i+1,k)=L0;
       dL(i+1,k)=L0_1;
    end %i
end %k
