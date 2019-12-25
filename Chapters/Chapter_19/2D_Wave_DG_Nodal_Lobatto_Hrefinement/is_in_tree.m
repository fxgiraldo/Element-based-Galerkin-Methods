%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [ptr,sfc,m] = is_in_tree(sfc,tm,tc,ptr,m)

if (tm(ptr)==1)
    m=m+1;
    sfc(m) = ptr;
else
    if(tc(1,ptr)~=0)
       for i=1:4
          [ptr1,sfc,m] = is_in_tree(sfc,tm,tc,tc(i,ptr),m);
       end
    end
end
