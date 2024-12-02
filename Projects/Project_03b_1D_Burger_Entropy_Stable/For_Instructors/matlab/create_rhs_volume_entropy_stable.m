%---------------------------------------------------------------------%
%This function computes the Volume Integrals for Energy-Stable flux.
%Written by F.X. Giraldo on 5/2021
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_volume_entropy_stable(qp,intma,coord,npoin,nelem,ngl,wgl,dpsi,flux_method,cgdg_method)

if strcmp(cgdg_method,'strong')
    rhs = create_rhs_volume_strong_entropy_stable(qp,intma,coord,npoin,nelem,ngl,wgl,dpsi,flux_method);
elseif strcmp(cgdg_method,'weak')
    rhs = create_rhs_volume_weak_entropy_stable(qp,intma,coord,npoin,nelem,ngl,wgl,dpsi,flux_method);
end