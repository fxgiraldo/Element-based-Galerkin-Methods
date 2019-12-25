%---------------------------------------------------------------------%
%This function computes the 2D Mass Matrix on Quadrilaterals using 
%tensor product of 1D basis functions and Inexact Integration.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [Mmatrix, Mmatrix_global] = create_Mmatrix2d_dg_TensorProduct_inexact_all(jac,intma,iperiodic,npoinT,nelemT,ngl,face,nface)

%Initialize
Mmatrix=zeros(5*nelemT,ngl,ngl);
Mmatrix_global=zeros(npoinT,1);
% 
% for is=1:nsfc
%     e = sfc(is);
for e=1:nelemT
   %Do LGL Integration
    for j=1:ngl
        for i=1:ngl 
            ip=iperiodic(intma(e,i,j));
            wq=jac(e,i,j);
            Mmatrix(e,i,j)=Mmatrix(e,i,j) + wq;
            Mmatrix_global(ip)=Mmatrix_global(ip) + wq;
        end %i
    end %j
end %e

% %gather from hanging nodes
% Mmatrix_global = gather_from_hanging_nodes_CG_DG(Mmatrix,Mmatrix_global,intma,face,ffc,nffc,ngl);
Mmatrix_global = gather_from_hanging_nodes_CG_DG_all(Mmatrix,Mmatrix_global,intma,face,nface,ngl);

%Periodicity
for i=1:npoinT
 j=iperiodic(i);
    if (i ~= j) 
        Mmatrix_global(i)=1;
    end
end %i 
      