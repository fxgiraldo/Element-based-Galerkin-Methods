%---------------------------------------------------------------------%
%This function computes the 2D Mass Matrix on Quadrilaterals using 
%tensor product of 1D basis functions and Inexact Integration.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [Mmatrix] = create_Mmatrix2d_TensorProduct_inexact(jac,intma,nelem,ngl)

%Initialize
Mmatrix=zeros(nelem,ngl,ngl);
Mmatrix=jac;

% for e=1:nelem
%    %Do LGL Integration
%     for j=1:ngl
%         for i=1:ngl 
%             wq=jac(e,i,j);
%             Mmatrix(e,i,j)=wq;
%         end %i
%     end %j
% end %e
      