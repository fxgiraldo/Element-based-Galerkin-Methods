%---------------------------------------------------------------------%
%This code computes the Legendre-Gauss-Lobatto points and weights
%which are the roots of the Lobatto Polynomials.
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [xgl,wgl] = legendre_gauss_lobatto(P)

p=P-1; %Order of the Polynomials
ph=floor( (p+1)/2 );

for i=1:ph
   x=cos( (2*i-1)*pi/(2*p+1) );
   for k=1:20
      [L0,L0_1,L0_2]=legendre_poly(p,x); %Compute Nth order Derivatives of Legendre Polys
      
      %Get new Newton Iteration
      dx=-(1-x^2)*L0_1/(-2*x*L0_1 + (1-x^2)*L0_2);
      x=x+dx;
      if (abs(dx) < 1.0e-20) 
         break
      end
   end
   xgl(p+2-i)=x;
   wgl(p+2-i)=2/(p*(p+1)*L0^2);
end

%Check for Zero Root
if (p+1 ~= 2*ph)
   x=0;
   [L0,L0_1,L0_2]=legendre_poly(p,x);
   xgl(ph+1)=x;
   wgl(ph+1)=2/(p*(p+1)*L0^2);
end
   
%Find remainder of roots via symmetry
for i=1:ph
   xgl(i)=-xgl(p+2-i);
   wgl(i)=+wgl(p+2-i);
end
   


      
