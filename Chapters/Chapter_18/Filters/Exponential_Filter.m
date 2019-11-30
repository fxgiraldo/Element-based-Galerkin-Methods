%---------------------------------------------------------------------%
%This code computes the Exponential Filter Function
%Written by F.X. Giraldo on 1/2016
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [sigma] = Exponential_Filter(N,alpha,F,k,s,kmax)

%Constants
sigma=zeros(kmax,1);

%Compute Weight
for i=1:kmax
    sigma(i)=exp(  - alpha*( (k(i))/N )^F );
end



      
