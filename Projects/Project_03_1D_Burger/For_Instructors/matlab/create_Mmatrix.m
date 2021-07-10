%---------------------------------------------------------------------%
%This function computes the global mass matrix for either CG or DG.
%Written by F.X. Giraldo on 8/2015
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Mmatrix = create_Mmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi)

%Initialize
Mmatrix=zeros(npoin,npoin);

for e=1:nelem

    %Store Coordinates
    for i=1:ngl
        x(i)=coord(i,e);
    end

    dx=x(ngl)-x(1);
    jac=dx/2;

    %Do LGL Integration
    for l=1:nq
        wq=wnq(l)*jac;
        for i=1:ngl
            I=intma(i,e);
            h_i=psi(i,l);
            for j=1:ngl
                J=intma(j,e);
                h_j=psi(j,l);
                Mmatrix(I,J)=Mmatrix(I,J) + wq*h_i*h_j;
            end %j
        end %i
    end %l
end %e



      