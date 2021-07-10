%---------------------------------------------------------------------%
%This function computes the RK Coefficients.
%Written by F.X. Giraldo on May 7, 2012
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [a0,a1,beta] = compute_ti_coefficients(kstages)

a0=zeros(kstages);
a1=zeros(kstages);
beta=zeros(kstages);

if kstages == 1 %RK1-SSP
    a0(1)=1; a1(1)=0; beta(1)=1;
elseif kstages == 2 %RK2-SSP
    a0(1)=1; a1(1)=0; beta(1)=1; %SSP
    a0(2)=1.0/2.0; a1(2)=1.0/2.0; beta(2)=1.0/2.0;
%     a0(1)=1; a1(1)=0; beta(1)=0.5; %Non-SSP
%     a0(1)=1; a1(1)=0; beta(1)=1;
elseif kstages == 3 %RK3-SSP
    a0(1)=1; a1(1)=0; beta(1)=1; 
    a0(2)=3.0/4.0; a1(2)=1.0/4.0; beta(2)=1.0/4.0;
    a0(3)=1.0/3.0; a1(3)=2.0/3.0; beta(3)=2.0/3.0;
elseif kstages == 4 %RK4-Non-SSP
    a0(1)=1; a1(1)=0; beta(1)=1.0/2.0;
    a0(2)=1; a1(2)=0; beta(2)=1.0/2.0;
    a0(3)=1; a1(3)=0; beta(3)=1.0;
    a0(4)=1; a1(4)=0; beta(4)=0;
end