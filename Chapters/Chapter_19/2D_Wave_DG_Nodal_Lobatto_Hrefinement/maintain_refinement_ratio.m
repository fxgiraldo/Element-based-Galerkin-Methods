%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [nref,iref,ncor,icor] = maintain_refinement_ratio(nref,iref,ncor,icor,jeside,face)
%UNFINISHED
for i=1:nref
    parent = iref(i);
    [cr,pref,npref] = check_ratio(parent,jeside,face);
    temp = iref;
    iref = npref;
    for j=1:nref
        iref
    end
end


end
