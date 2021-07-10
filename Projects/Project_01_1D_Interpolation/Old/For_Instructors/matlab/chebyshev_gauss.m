%---------------------------------------------------------------------%
%This code computes the Chebyshev points and weights
%Written by F.X. Giraldo on 7/2014
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [xgl,wgl] = chebyshev_gauss(P)

p=P-1; %Order of the Polynomials

for i=1:P
   xgl(i)=cos( (2*i-1)*pi/(2*P) ); 
   wgl(i)=pi/P;
end

   


      
