%---------------------------------------------------------------------%
%This function unites DPSI from the left and right portions of the domain
%in order to build one DPSI for each element.
%Written by F.X. Giraldo on January 19, 2024.
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function dpsi=unite_dpsi(dpsil,dpsir,ngl,nelem,nelem_LGL)

%Store Basis Function derivatives in one element-wise array
dpsi=zeros(ngl(nelem),ngl(nelem),nelem);
for e=1:nelem_LGL
    for j=1:ngl(e)
        for i=1:ngl(e)
            dpsi(i,j,e)=dpsil(i,j);
        end
    end
end
for e=nelem_LGL+1:nelem
    for j=1:ngl(e)
        for i=1:ngl(e)
            dpsi(i,j,e)=dpsir(i,j);
        end
    end
end