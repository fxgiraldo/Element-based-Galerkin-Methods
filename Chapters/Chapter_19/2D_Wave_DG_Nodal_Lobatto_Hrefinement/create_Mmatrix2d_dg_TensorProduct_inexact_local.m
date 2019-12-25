%---------------------------------------------------------------------%
%This function computes the 2D Mass Matrix on Quadrilaterals using 
%tensor product of 1D basis functions and Inexact Integration.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [Mmatrix,Mmatrix_inv] = create_Mmatrix2d_dg_TensorProduct_inexact_local(Mmatrix,Mmatrix_inv,jac,ngl,new_el,nnew)
 
for is=1:nnew
    e=new_el(is);
   %Do LGL Integration
    for j=1:ngl
        for i=1:ngl 
            wq=jac(e,i,j);
            Mmatrix(e,i,j)=Mmatrix(e,i,j) + wq;
        end %i
    end %j
    Mmatrix_inv(e,:,:)=1/Mmatrix(e,:,:);
end %e


end
      