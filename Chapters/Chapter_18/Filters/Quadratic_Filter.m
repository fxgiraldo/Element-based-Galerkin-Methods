%---------------------------------------------------------------------%
%This code computes the Quadratic Filter Function
%Written by F.X. Giraldo on 1/2016
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [sigma] = Quadratic_Filter(N,alpha,F,k,s,kmax)

%Constants
sigma=zeros(kmax,1);

%Compute Weight
for i=1:kmax
    if (k(i) <= s)
        sigma(i)=1;
    elseif (k(i) > s)
        x=(k(i)-s)/(N-s);
        sigma(i)=1 - alpha*x^F;
    end
end



      
