%----------------------------------------------------------------------%
%This subroutine builds the RHS vector for the Weak Form CG/DG Methods
%on Quadrilateral Elements for the 2D Shallow Water Equations.
%Written by Francis X. Giraldo on 11/2022
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5502
%----------------------------------------------------------------------%
function rhs = create_rhs(q,ksi_x,ksi_y,eta_x,eta_y,jac,...
		 wnq,psi,dpsi,intma,iperiodic,Mmatrix,face,normals,jac_face,...
         mapL,mapR,npoin,nelem,nface,ngl,space_method,flux_method,form_method)

%------------Students Add Your Routines Here---------------%
%------------Students Add Your Routines Here---------------%
%Compute Volume Integral: Differentiation matrix contribution
rhs = create_rhs_volume(q,ksi_x,ksi_y,eta_x,eta_y,jac,...
      wnq,psi,dpsi,intma,iperiodic,npoin,nelem,ngl,form_method);
%Construct Communicator: DSS for CG / Fluxes for DG
rhs = create_rhs_flux(rhs,q,face,normals,jac_face, ...
      wnq,nface,ngl,mapL,mapR,intma,space_method,flux_method,form_method);
%------------Students Add Your Routines Here---------------%
%------------Students Add Your Routines Here---------------%

%Multiply by Inverse Mass
for i=1:npoin
for l=1:3
    rhs(l,i) = rhs(l,i)/Mmatrix(i);
end
end

%Periodicity
for i=1:npoin
    for l=1:3
        rhs(l,i)=rhs(l,iperiodic(i));
    end
end
