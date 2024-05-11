%---------------------------------------------------------------------%
%This function computes the Mass and Energy for Burgers Equation
%Written by F.X. Giraldo on 11/14/2018 from ECMWF
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [m1] = compute_Mass_and_Energy(qp,intma,jac,wnq,nelem,ngl)

%Compute Mass and Energy
m1=0;
for e=1:nelem
    for l=1:ngl(e)
      wq=wnq(l,e)*jac(e);
      for j=1:ngl(e)
        jp=intma(j,e);
        h=(qp(jp));
        m1=m1 + wq*h;
      end
    end
end