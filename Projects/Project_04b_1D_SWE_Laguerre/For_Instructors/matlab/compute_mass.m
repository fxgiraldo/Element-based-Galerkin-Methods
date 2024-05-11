%---------------------------------------------------------------------%
%This function computes the mass and energy.
%Written by F.X. Giraldo on January 19, 2024.
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [mass,energy]=compute_mass(qp,qb,intma,ngl,nelem,wnq,jac)

%Compute Initial Mass
m1=0; m2=0;
for e=1:nelem
    for i=1:ngl(e)
        wq=wnq(i,e)*jac(e);
        I=intma(i,e);
        h=qp(I,1)+qb(I);
        U=qp(I,2);
        m1=m1 + wq*h;
        m2=m2 + wq*U;
    end
end
mass=m1;
energy=m2;
