%---------------------------------------------------------------------%
%This function computes the 2D Laplacian Matrix on Quadrilaterals.
%Written by F.X. Giraldo on May 29, 2024
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [Mmatrix,Lmatrix] = create_Global_Matrices(intma,jac,wnq,ksi_x,ksi_y,eta_x,...         
                   eta_y,psi,dpsi,iperiodic,npoin,nelem,ngl,nq)

if (ngl==nq)
    [Mmatrix,Lmatrix] = create_Global_Matrices_Inexact(intma,jac,wnq,ksi_x,ksi_y,eta_x,...         
                   eta_y,psi,dpsi,iperiodic,npoin,nelem,ngl,nq);
elseif (ngl~=nq)
    [Mmatrix,Lmatrix] = create_Global_Matrices_Exact(intma,jac,wnq,ksi_x,ksi_y,eta_x,...         
                   eta_y,psi,dpsi,iperiodic,npoin,nelem,ngl,nq);
end

      
