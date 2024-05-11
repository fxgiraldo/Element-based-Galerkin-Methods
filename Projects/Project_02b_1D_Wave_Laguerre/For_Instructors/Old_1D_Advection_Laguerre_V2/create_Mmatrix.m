%---------------------------------------------------------------------%
%This function computes the global mass matrix for either CG or DG.
%Written by F.X. Giraldo on 8/2015
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Mmatrix = create_Mmatrix(intma,jac,wnq,npoin,nelem,ngl)

%Initialize
Mmatrix=zeros(npoin,1);

for e=1:nelem
    for i=1:ngl(e)
        wq=wnq(i,e)*jac(e);
        I=intma(i,e);
        Mmatrix(I)=Mmatrix(I) + wq;
    end %i
end %e



      