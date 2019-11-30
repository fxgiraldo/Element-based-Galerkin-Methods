%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [M,D,F] = global_matrices(Me,De,Fe,intma,coord,npoin,nelem,ngl,periodicity)

%Initialize
M=zeros(npoin,npoin);
D=zeros(npoin,npoin);
F=zeros(npoin,npoin);

%Form Global Mass and Differentiation Matrices
for e=1:nelem
    for i=1:ngl
        I=intma(i,e);
        x(i)=coord(I);
        inode(i)=periodicity(I);
    end
    jac=0.5*(x(ngl)-x(1));
    for j=1:ngl
        J=inode(j);
        for i=1:ngl
            I=inode(i);
            M(I,J)=M(I,J) + jac*Me(i,j);
            D(I,J)=D(I,J) + De(i,j);
            F(I,J)=F(I,J) + Fe(i,j);
        end
    end
end 

if periodicity(npoin) == periodicity(1)
    M(npoin,npoin)=1;
end
