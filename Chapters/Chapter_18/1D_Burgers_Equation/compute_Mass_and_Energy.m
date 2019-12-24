%---------------------------------------------------------------------%
%This function computes the Mass and Energy for Burgers Equation
%Written by F.X. Giraldo on 11/14/2018 from ECMWF
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [m1,m2] = compute_Mass_and_Energy(qp,intma,coord,wnq,psi,nelem,ngl,nq)

%Compute Mass and Energy
m1=0; m2=0;
for e=1:nelem
   dx=coord(ngl,e)-coord(1,e);
   jac=dx/2;
    for l=1:nq
      wq=wnq(l)*jac;
      for j=1:ngl
         jp=intma(j,e);
         h=(qp(jp,1));
         U=0.5*h^2;
         m1=m1 + wq*(h)*psi(j,l);
         m2=m2 + wq*(U)*psi(j,l);
      end
    end
end