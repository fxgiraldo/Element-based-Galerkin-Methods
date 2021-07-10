%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on April 22, 2021
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [M,D] = Matrix_DSS(Me,De,u,intma,periodicity,ngl,nelem,npoin)

%Form Global Matrices
M=zeros(npoin,npoin);
D=zeros(npoin,npoin);

for e=1:nelem
    for i=1:ngl
        ip=periodicity(intma(i,e));
        for j=1:ngl
            jp=periodicity(intma(j,e));
            M(ip,jp)=M(ip,jp) + Me(i,j,e);
            D(ip,jp)=D(ip,jp) + u*De(i,j);
        end
    end
end


      