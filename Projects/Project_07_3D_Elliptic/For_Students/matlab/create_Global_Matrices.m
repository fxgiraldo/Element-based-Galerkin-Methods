%---------------------------------------------------------------------%
%This function computes the 3D Laplacian Matrix on Quadrilaterals.
%Written by F.X. Giraldo on June, 2024
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [Mmatrix,Lmatrix] = create_Global_Matrices(intma,jac,wnq,ksi_x,ksi_y,ksi_z, ...
                             eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,psi,dpsi, ...
                             npoin,nelem,ngl,nq)

if (ngl==nq)
    [Mmatrix,Lmatrix] = create_Global_Matrices_Inexact(intma,jac,wnq,ksi_x,ksi_y,ksi_z, ...
                       eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,psi,dpsi, ...
                       npoin,nelem,ngl,nq);

elseif (ngl~=nq)
    [Mmatrix,Lmatrix] = create_Global_Matrices_Exact(intma,jac,wnq,ksi_x,ksi_y,ksi_z, ...
                       eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,psi,dpsi, ...
                       npoin,nelem,ngl,nq);

end

      
