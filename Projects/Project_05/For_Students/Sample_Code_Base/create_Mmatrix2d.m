%---------------------------------------------------------------------%
%This function computes the 2D Mass Matrix on Quadrilaterals using 
%tensor product of 1D basis functions and Inexact Integration.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Mmatrix = create_Mmatrix2d(jac,wnq,intma,iperiodic,npoin,nelem,ngl)

%Initialize
Mmatrix=zeros(npoin,1);

for e=1:nelem
   %Do LGL Integration
    for j=1:ngl
        for i=1:ngl 
            I=iperiodic(intma(e,i,j));
            wq=wnq(i)*wnq(j)*jac(e,i,j);
            Mmatrix(I)=Mmatrix(I) + wq;
        end %i
    end %j
end %e

%Periodicity
for i=1:npoin
 j=iperiodic(i);
    if (i ~= j) 
        Mmatrix(i)=1;
    end
end %i 
  

      