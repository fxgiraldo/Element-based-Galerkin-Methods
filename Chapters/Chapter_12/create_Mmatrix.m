%---------------------------------------------------------------------%
%This function computes the 2D Mass Matrix on Quadrilaterals using 
%tensor product of 1D basis functions and Inexact Integration.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Mmatrix = create_Mmatrix(jac,intma,psi,npoin,nelem,ngl,nq)

%Initialize
Mmatrix=zeros(npoin,npoin);

for e=1:nelem
   
   %Do LGL Integration
   for k2=1:nq
   for k1=1:nq
        wq=jac(e,k1,k2);
        for j2=1:ngl
        for j1=1:ngl
            jp=intma(e,j1,j2);
            h_j=psi(j1,k1)*psi(j2,k2);
            for i2=1:ngl
            for i1=1:ngl 
                ip=intma(e,i1,i2);
                h_i=psi(i1,k1)*psi(i2,k2);
                Mmatrix(ip,jp)=Mmatrix(ip,jp) + wq*h_i*h_j;
            end
            end
        end
        end
   end 
   end 
end %e  

      