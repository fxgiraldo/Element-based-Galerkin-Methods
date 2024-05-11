function w = vandeven_modal(F,N,cutoff,M)
% w = vandeven_modal(F,N,cutoff)
% Construct a generalized Boyd-Vandeven filter



%Default Boyd-Vandeven filter
if (nargin == 2)
    cutoff = 2/3;
    M = N;
end


if (nargin == 3)
    M = N;
end


if ( (cutoff > 1) | (cutoff <= 0 ))
    error('cutoff must lie between 0 and 1')
end



Np1 = N + 1;
w = zeros(Np1,1);            %Filter weights
s = ceil(cutoff*M);
rootF = sqrt(F);

for k = 0:(M -1)
    k1 = k + 1; 
    w(k1) = 1;
end    


for k = s:M
   k1 = k + 1; 
   x = (k - s)/(N - s);
   Omega  = abs(x) - .5;
   xlog = log(1 - 4.* Omega.^2);
   diff = abs(x - .5);
   if (diff < 1e-10)
       square_root = 1;
   else
       square_root =  sqrt( - log(1 - 4.* Omega.^2)/(4.*Omega.^2));;
   end
   
   
   
  
   sigma = .5 * erfc ( 2*rootF .* Omega.* square_root);
   w(k1) = sigma;
end





