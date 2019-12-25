%---------------------------------------------------------------------%
%This code computes the Lagrange Polynomials at given x locations
%Written by M.A. Kopera on 10/2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%

function L=lagrange_poly(x,xgl)

nx = size(x,2);
ngl=size(xgl,2);
L=ones(ngl,nx);
   for i=1:ngl
      for j=1:ngl
         if (i~=j)
            L(i,:)=L(i,:).*(x-xgl(j))/(xgl(i)-xgl(j));
         end
      end
   end
end