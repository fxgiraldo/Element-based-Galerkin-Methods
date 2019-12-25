%---------------------------------------------------------------------%
%This subroutine creates the space filling curve out of the element tree
% note that nelem should be the original number of elements at the top
% level
% elmptr is the element pointer array which defines the element tree
% recursively
% nsfc is the length of sfc
%Written by M.A. Kopera on 10/2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%

function [sfc,nsfc] = space_filling_curve(nelem,elmptr)

clear sfc;
m=0;
for e=1:nelem
    if(elmptr(e)==0)
       m=m+1;
       sfc(m)=e;
    else
       el=elmptr(e);
       for k=1:4
        m=m+1;
        sfc(m)=el;
        el=el+1;
       end
    end
end
nsfc=m;

end