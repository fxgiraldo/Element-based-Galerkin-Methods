%---------------------------------------------------------------------%
%This function computes the RHS Volume terms.
%Written by F.X. Giraldo on 1/2024
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_volume(qp,qb,intma,jac,npoin,nelem,ngl,wnq,dpsi,gravity,delta_nl,Igamma,icase,form_method)

if strcmp(form_method,'strong')
    rhs = create_rhs_volume_strong(qp,qb,intma,jac,npoin,nelem,ngl,wnq,dpsi,gravity,delta_nl,Igamma,icase);
elseif strcmp(form_method,'weak')
    rhs = create_rhs_volume_weak(qp,qb,intma,jac,npoin,nelem,ngl,wnq,dpsi,gravity,delta_nl,Igamma,icase);
end
