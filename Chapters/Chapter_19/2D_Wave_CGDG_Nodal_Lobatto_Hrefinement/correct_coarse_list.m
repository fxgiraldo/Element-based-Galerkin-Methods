%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [icorc,ncorc] = correct_coarse_list(icor,ncor,tp,tc)

%sort
icor = sort(icor);
icorc = 0;
ncorc = 0;
i=1;
m=1;
k=0;
while i<=ncor

    if icor(i)==tc(m,tp(icor(i))) 
        i=i+1;
        m=m+1;
    else

        m=1;
        i=i+1;
        continue
    end
    
    if m==5
        k=k+1;
        icorc(k) = icor(i-1);
        m=1;
    end
end

ncorc = k;
end
