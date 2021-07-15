%----------------------------------------------------------------------%
%This subroutine builds the RHS vector for the Weak Form CG/DG Methods
%on Quadrilateral Elements for the 2D Advection Equation.
%Written by Francis X. Giraldo on 5/2021
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5502
%----------------------------------------------------------------------%
function rhs = construct_RHS_vector(q,u,v,ksi_x,ksi_y,eta_x,eta_y,jac,...
		 wnq,dpsi,intma,iperiodic,Mmatrix,face,normals,jac_face,...
         mapL,mapR,npoin,nelem,nface,ngl,space_method,flux_method)

%------------Students Add Your Routines Here---------------%
%Compute Volume Integral: Differentiation matrix contribution
rhs = compute_rhs_differentiation(q,u,v,ksi_x,ksi_y,eta_x,eta_y,jac,...
    wnq,dpsi,intma,iperiodic,npoin,nelem,ngl);
%Construct Communicator: DSS for CG / Fluxes for DG
if strcmp(space_method,'dg')
    rhs = compute_rhs_flux(rhs,q,u,v,face,normals,jac_face, ...
          wnq,nface,ngl,mapL,mapR,intma,flux_method);
end %if
%------------Students Add Your Routines Here---------------%

%Multiply by Inverse Mass
rhs = rhs./Mmatrix;

%Periodicity
for i=1:npoin
    rhs(i)=rhs(iperiodic(i));
end
