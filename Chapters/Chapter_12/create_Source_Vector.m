%---------------------------------------------------------------------%
%This function computes the 2D Mass Matrix on Quadrilaterals using 
%tensor product of 1D basis functions and Inexact Integration.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function R = create_Source_Vector(jac,intma,psi,npoin,nelem,ngl,nq,fe)

%Initialize
R=zeros(npoin,1);

for e=1:nelem
   
   %Do LGL Integration
   for k2=1:nq
   for k1=1:nq
        wq=jac(e,k1,k2);
        
        %Interpolate FE to Quadrature Points
        f_k=0;
        for j2=1:ngl
        for j1=1:ngl
            jp=intma(e,j1,j2);
            h_j=psi(j1,k1)*psi(j2,k2);
            f_k=f_k + h_j*fe(jp);
        end
        end
        
        for i2=1:ngl
        for i1=1:ngl 
            ip=intma(e,i1,i2);
            h_i=psi(i1,k1)*psi(i2,k2);
            R(ip)=R(ip) + wq*h_i*f_k;
        end
        end
   end 
   end 
end %e  

      