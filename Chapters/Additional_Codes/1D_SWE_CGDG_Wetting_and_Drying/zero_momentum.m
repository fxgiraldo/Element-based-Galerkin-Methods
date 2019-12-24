%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Maths
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function qp = zero_momentum(qp,qb,nelem,ngl,h_eps)

%Integrate Divergence of Flux (Weak Form)
for e=1:nelem
   
   %Store Coordinates
   for i=1:ngl
      p_k=qp(1,i,e) + qb(i,e);
      if (p_k <= h_eps)
          qp(2,e)=0;
      end
   end
end %e
